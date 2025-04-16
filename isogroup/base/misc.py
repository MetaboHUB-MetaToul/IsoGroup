import re
from isocor.base import LabelledChemical

class Misc:

    @staticmethod
    def _parse_strtracer(str_tracer):
        """Parse the tracer code.

        Args:
            str_tracer (str): tracer code (e.g. "13C")

        Returns:
            tuple
                - (str) tracer element (e.g. "C")
                - (int) tracer index in :py:attr:`~data_isotopes`
        """
        try:
            tracer = re.search(r'(\d*)([A-Z][a-z]*)', str_tracer)
            count = int(tracer.group(1))
            tracer_el = tracer.group(2)
        except (ValueError, AttributeError):
            raise ValueError("Invalid tracer code: '{}'."
                             " Please check your inputs.".format(str_tracer))
        best_diff = float("inf")
        idx_tracer = None
        unexpected_msg = "Unexpected tracer code. Are you sure this isotope is "\
                         "in data_isotopes? '{}'".format(str_tracer)
        assert tracer_el in LabelledChemical.DEFAULT_ISODATA, unexpected_msg
        for i, mass in enumerate(LabelledChemical.DEFAULT_ISODATA[tracer_el]["mass"]):
            test_diff = abs(mass - count)
            if test_diff < best_diff:
                best_diff = test_diff
                idx_tracer = i
        assert best_diff < 0.5, unexpected_msg
        return (tracer_el, idx_tracer)