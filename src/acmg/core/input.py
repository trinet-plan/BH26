"""Lossless demo audit alongside a strict, allowlisted evaluation input."""

import hashlib
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree

from acmg.core.models import Variant


IDENTITY_FIELDS = ("GENE", "TRANSCRIPT", "HGVSC", "HGVSP", "CLNVARIATIONID")
ANNOTATION_FIELDS = (
    "GNOMAD_AF", "GNOMAD_AF_EAS", "GNOMAD_AF_LATINO", "SG10K_AF", "TOMMO_AF", "HGVD_AF",
    "AM_PATHOGENICITY", "AM_CLASS", "AG_SPLICING_SCORE", "AG_IMPACT_RAW", "AG_IMPACT_PHRED",
)
XLSX_ANNOTATION_FIELDS = (
    "am_pathogenicity", "am_class", "am_protein_variant", "am_transcript",
    "ag_splicing_score", "ag_impact_raw_score", "ag_impact_PHRED",
)
_CELL_REF = re.compile(r"([A-Z]+)")
_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_info(value):
    result = {}
    for item in value.split(";"):
        if item and item != ".":
            key, sep, val = item.partition("=")
            if key in result:
                raise ValueError(f"Duplicate INFO key: {key}")
            result[key] = val if sep else True
    return result


def _column_index(reference):
    letters = _CELL_REF.match(reference).group(1)
    result = 0
    for letter in letters:
        result = result * 26 + ord(letter) - 64
    return result - 1


def read_xlsx_rows(path):
    """Read the first worksheet as typed rows using only the Python standard library."""
    path = Path(path)
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", _NS):
                shared.append("".join(node.text or "" for node in item.iterfind(".//m:t", _NS)))
        root = ElementTree.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    rows = []
    for row in root.findall(".//m:sheetData/m:row", _NS):
        values = {}
        for cell in row.findall("m:c", _NS):
            index = _column_index(cell.attrib["r"])
            cell_type = cell.attrib.get("t")
            value_node = cell.find("m:v", _NS)
            if cell_type == "inlineStr":
                value = "".join(node.text or "" for node in cell.iterfind(".//m:t", _NS))
            elif value_node is None:
                value = None
            elif cell_type == "s":
                value = shared[int(value_node.text)]
            elif cell_type == "b":
                value = value_node.text == "1"
            else:
                number = value_node.text
                value = float(number) if "." in number or "e" in number.lower() else int(number)
            values[index] = value
        width = max(values, default=-1) + 1
        rows.append([values.get(index) for index in range(width)])
    if not rows:
        return []
    headers = rows[0]
    while headers and headers[-1] is None:
        headers.pop()
    if any(not isinstance(value, str) or not value for value in headers):
        raise ValueError(f"Spreadsheet has a blank or invalid header: {path}")
    result = []
    for row_number, values in enumerate(rows[1:], 2):
        padded = (values + [None] * len(headers))[:len(headers)]
        result.append({"row": row_number, **dict(zip(headers, padded, strict=True))})
    return result


def attach_demo_spreadsheet(records, path):
    """Audit spreadsheet annotations without promoting labels into decision evidence."""
    path = Path(path)
    digest = sha256_file(path)
    rows = read_xlsx_rows(path)
    by_key = {}
    for row in rows:
        key = (row.get("case"), row.get("variant_id"))
        if not all(key) or key in by_key:
            raise ValueError(f"Invalid or duplicate spreadsheet key at row {row['row']}")
        by_key[key] = row
    for record in records:
        key = (record["source"]["case_id"], record["source"]["variant_id"])
        row = by_key.pop(key, None)
        if row is None:
            record["issues"].append("MISSING_SPREADSHEET_ANNOTATION")
            continue
        record["spreadsheet_source"] = {
            "file": path.name, "sha256": digest, "sheet": "Sheet1", "row": row["row"]
        }
        record["spreadsheet_annotations"] = {
            key: row.get(key) for key in XLSX_ANNOTATION_FIELDS if row.get(key) is not None
        }
        comparisons = {
            "gene": record["identity"].get("GENE"),
            "hgvsc": record["identity"].get("HGVSC"),
            "hgvsp": record["identity"].get("HGVSP"),
            "chrom": int(record["raw_variant"]["chrom"]),
            "pos": int(record["raw_variant"]["pos"]),
            "vcf_ref": record["raw_variant"]["ref"],
            "vcf_alt": record["raw_variant"]["alt"],
        }
        if any(row.get(field) != value for field, value in comparisons.items()):
            record["issues"].append("VCF_SPREADSHEET_IDENTITY_MISMATCH")
    if by_key:
        raise ValueError(f"Spreadsheet contains records absent from VCFs: {sorted(by_key)}")
    return records


def audit_vcf(path, case_id=None, assembly=None):
    """Retain every source row and every ALT, even malformed or unsupported records."""
    path = Path(path)
    digest = sha256_file(path)
    records = []
    reference = None
    caveats = []
    header_seen = False
    for line_no, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if line.startswith("##reference="):
            reference = line.partition("=")[2]
        if line.startswith("##caveat="):
            caveats.append(line.partition("=")[2])
        if line.startswith("#CHROM\t"):
            header_seen = True
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        errors = []
        if not header_seen:
            errors.append("MISSING_VCF_COLUMN_HEADER")
        if len(fields) < 8:
            errors.append("INVALID_COLUMN_COUNT")
        info = {}
        if len(fields) >= 8:
            try:
                info = parse_info(fields[7])
            except ValueError as exc:
                errors.append(str(exc))
        alts = fields[4].split(",") if len(fields) >= 5 else [None]
        for alt_index, alt in enumerate(alts, 1):
            issues = errors.copy()
            source = {
                "file": path.name, "sha256": digest, "line": line_no,
                "case_id": case_id or path.stem, "variant_id": fields[2] if len(fields) > 2 else None,
                "alt_index": alt_index,
            }
            raw_variant = {
                "assembly": assembly or reference,
                "chrom": fields[0], "pos": fields[1] if len(fields) > 1 else None,
                "ref": fields[3] if len(fields) > 3 else None, "alt": alt,
            }
            parsed_variant = None
            if assembly and reference and assembly != reference:
                issues.append("ASSEMBLY_CONFLICT")
            try:
                parsed_variant = Variant(**{**raw_variant, "pos": int(raw_variant["pos"])}).to_dict()
            except (ValueError, TypeError) as exc:
                issues.append(f"INVALID_VARIANT: {exc}")
            identity = {key: info[key] for key in IDENTITY_FIELDS if key in info}
            if str(identity.get("HGVSC", "")).find("+") >= 0 and "HGVSP" in identity:
                if identity["HGVSP"] not in {"p.?", "."}:
                    issues.append("CHECK_INTRONIC_HGVSC_VS_PROTEIN_ANNOTATION")
            records.append({
                "record_id": f"{source['case_id']}:{line_no}:{alt_index}",
                "source": source, "raw_variant": raw_variant, "identity": identity,
                "parsed_variant": parsed_variant,
                "identity_status": "PENDING", "issues": issues,
                "source_caveats": caveats.copy(),
                "original_annotations": {key: info[key] for key in ANNOTATION_FIELDS if key in info},
                "source_record": line,
                "resolution": None,
            })
    return records


def audit_demo(directory):
    directory = Path(directory)
    records = []
    for case in range(1, 5):
        path = directory / f"case{case}_variants_v2.vcf"
        if not path.is_file():
            raise ValueError(f"Missing demo input: {path}")
        records.extend(audit_vcf(path, case_id=f"case{case}"))
    # The same transcript HGVS appearing with different genomic alleles is an identity conflict.
    by_hgvs = {}
    for record in records:
        identity = record["identity"]
        key = (identity.get("TRANSCRIPT"), identity.get("HGVSC"))
        if all(key) and record["parsed_variant"]:
            by_hgvs.setdefault(key, set()).add(Variant(**record["parsed_variant"]).key)
    for record in records:
        identity = record["identity"]
        key = (identity.get("TRANSCRIPT"), identity.get("HGVSC"))
        if len(by_hgvs.get(key, set())) > 1:
            record["issues"].append("CONFLICTING_GENOMIC_REPRESENTATIONS_FOR_HGVS")
    spreadsheet = directory / "annotation_alphamissense_alphagenome_v1.xlsx"
    if not spreadsheet.is_file():
        raise ValueError(f"Missing demo input: {spreadsheet}")
    return attach_demo_spreadsheet(records, spreadsheet)
