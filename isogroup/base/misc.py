from __future__ import annotations
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
    

    @staticmethod
    def get_atomic_mass(element: str) -> float|None:
        """
        Returns the atomic mass of the given element.
        Args:
            element (str): Chemical element symbol (e.g. "C", "H", "N", "O").
        Returns:
            float|None: Atomic mass of the element, or None if the element is not found
        """
        if element in LabelledChemical.DEFAULT_ISODATA:
            return LabelledChemical.DEFAULT_ISODATA[element]["mass"][0]
        return None
    

    @staticmethod
    def calculate_mzshift(tracer: str) -> float:
        """
        Calculate the m/z shift for a given tracer (e.g. "13C").
        Args:
            tracer (str): Tracer code (e.g. "13C").
        Returns:
            float: m/z shift corresponding to the tracer.
        """
        tracer_element, tracer_idx = Misc._parse_strtracer(tracer)

        # Mass of the tracer isotope
        tracer_mass = LabelledChemical.DEFAULT_ISODATA[tracer_element]["mass"][tracer_idx]

        # Mass of the most abundant natural isotope
        natural_mass = LabelledChemical.DEFAULT_ISODATA[tracer_element]["mass"][0]

        mz_shift = tracer_mass - natural_mass
        return mz_shift
    
    @staticmethod
    def get_max_isotopologues_for_mz(mz: float, tracer_element: str) -> int:
        """
        Returns the maximum number of isotopologues to consider based on the m/z value.
        This is a placeholder function and should be replaced with actual logic as needed.
        
        Args:
            mz (float): Mass-to-charge ratio of the feature.
            tracer_element (str): Tracer element symbol (e.g. "C", "N").
        
        Returns:
            int: Maximum number of isotopologues to consider.
        """
        element_mass = float(Misc.get_atomic_mass(tracer_element))
        if element_mass is None:
            raise ValueError(f"Unknown tracer element: {tracer_element}")
        if tracer_element == "C":
            factor = 0.7 # Approximation, empiric fraction based on the Seven Golden Rules
        elif tracer_element == "N":
            factor = 0.2
        elif tracer_element == "O":
            factor = 0.3
        else:
            raise NotImplementedError(f"Tracer {tracer_element} not implemented yet.")
        return max(1, int(factor * (mz / element_mass)))

    def calculate_isotopologue_index(candidate_mz, base_mz, mzshift_tracer):
        """
        Calculate the theoretical isotopologue index based on m/z values.
        Args:
            candidate_mz (float): m/z of the candidate isotopologue.
            base_mz (float): m/z of the base (unlabeled) feature.
            mzshift_tracer (float): m/z shift corresponding to the tracer.
        """
        return round((candidate_mz - base_mz) / mzshift_tracer)