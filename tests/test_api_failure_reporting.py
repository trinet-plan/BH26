"""A failed job must save why it failed.

The first real API run of MYBPC3 c.278delA failed because the demo VCF spells the deletion as
ALT ".", which AutomatedVariant refuses. What reached the saved result was
"ExceptionGroup: unhandled errors in a TaskGroup (1 sub-exception)" - the fact that something
broke, and nothing about what. run_pipeline() runs its providers in an anyio task group, so
every provider error arrives wrapped.
"""

import unittest

from api.main import describe_failure


def group(*exceptions):
    return ExceptionGroup("unhandled errors in a TaskGroup", list(exceptions))


class DescribeFailureTests(unittest.TestCase):

    REAL = ValueError("REF/ALT must be nonempty ACGT alleles; missing/SV ALT unsupported")

    def test_a_plain_exception_is_reported_as_itself(self):
        self.assertEqual(describe_failure(self.REAL),
                         {"type": "ValueError", "message": str(self.REAL)})

    def test_the_cause_is_dug_out_of_a_task_group(self):
        error = describe_failure(group(self.REAL))
        self.assertEqual(error["type"], "ValueError")
        self.assertEqual(error["message"], str(self.REAL))

    def test_the_wrapper_is_still_named(self):
        """Losing the cause was the bug; hiding that it was wrapped would be a second one."""
        self.assertEqual(describe_failure(group(self.REAL))["raised_as"], "ExceptionGroup")

    def test_a_cause_nested_several_groups_deep_is_found(self):
        self.assertEqual(describe_failure(group(group(group(self.REAL))))["message"],
                         str(self.REAL))

    def test_every_cause_is_kept_when_more_than_one_failed(self):
        """Reporting only the first would quietly drop the others."""
        error = describe_failure(group(self.REAL, KeyError("GENE")))
        self.assertEqual(error["type"], "ValueError")
        self.assertEqual([item["type"] for item in error["also_failed"]], ["KeyError"])

    def test_an_object_carrying_an_empty_exceptions_attribute_is_a_leaf(self):
        """Python forbids an empty ExceptionGroup, but `.exceptions` is just an attribute."""
        class Odd(Exception):
            exceptions = []
        self.assertEqual(describe_failure(Odd("nothing nested")),
                         {"type": "Odd", "message": "nothing nested"})


if __name__ == "__main__":
    unittest.main()
