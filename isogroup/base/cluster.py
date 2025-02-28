from typing import List, Union, Iterator, Self
from isogroup.base.feature import Feature
import re

class Cluster:

    def __init__(self, features: List[Feature], cluster_id = None):
        self.features = features
        self.cluster_id = cluster_id

    def __repr__(self) -> str:
        return f"Cluster({self.features})"
    
    
    def __len__(self) -> int:
        return len(self.features)

    def __add__(self, other: Union[Feature, Self]) -> Self:
        """
        Extend the cluster with either an extra feature or another cluster.
        :param other: Feature or Cluster
        :return: self
        """
        pass

    def __iadd__(self, other: Union[Feature, Self]) -> Self:
        """
        Extend the cluster with either an extra feature or another cluster.
        :param other: Feature or Cluster
        :return: self
        """
        pass

    def __contains__(self, item):
        pass

    def __iter__(self) -> Iterator[Feature]:
        return iter(self.features)
    
    
    def get_cluster(metabolite: str, sample: str) -> Self:
        """
        Returns the cluster of the metabolite in the given sample
        """
        pass


    def _get_element_number(self, element: str) -> int:
        """
        Returns the number of element tracer in the formula
        """
        formula = self.features[0].formula
        if formula is None:
            raise ValueError("Impossible to determine completeness without chemical formula.")
        
        # Extract the number of element tracer "C" in the formula   -- ## To be changed to be adaptable to any element tracer        
        element_number = re.findall(rf"{element}(\d+)", formula)

        if element in formula and not element_number:
            return 1
        elif element not in formula:
            return 0
            raise ValueError(f"The chemical formula does not contain the element '{element}'.")
        else:
            return int(element_number[0])
        

    @property
    def lowest_rt(self) -> float:
        return min([f.rt for f in self.features])

    @property
    def highest_rt(self) -> float:
        return max([f.rt for f in self.features])

    @property
    def lowest_mz(self) -> float:
        return min([f.mz for f in self.features])

    @property
    def highest_mz(self) -> float:
        return max([f.mz for f in self.features])

    @property
    def metabolite(self) -> str:
        """
        Returns the metabolite of the annotated features in the cluster
        """
        return self.features[0].metabolite
    
    
    @property
    def isotopologues(self) -> List[int]:
        """
        Returns the isotopologues of the annotated features in the cluster
        """
        return [f.isotopologue for f in self.features]
 

    @property
    def is_complete(self) -> bool:
        """
        Returns True if the cluster is complete
        """
        # Return an error if the cluster is not annotated
        if self.metabolite is None:
            raise ValueError("The cluster is not annotated. Please annotate the cluster first.")
        
        element_number = self._get_element_number("C")
        
        return (self.__len__()) == element_number + 1
        
    
    @property
    def missing_isotopologue(self) -> List[int]:
        """
        Returns a list of missing isotopologues in the cluster

        --> Adapt for more than one annotation (metabolite) in the cluster / 
        """
        # Check if the cluster is complete
        if self.is_complete:
            return []
        
        element_number = self._get_element_number("C")
        expected_isotopologues = element_number + 1
        current_isotopologues = set([item for sublist in self.isotopologues for item in sublist])

        missing_isotopologues = [i for i in range(expected_isotopologues) if i not in current_isotopologues]

        return missing_isotopologues
        
    
    @property
    def is_adduct(self) -> tuple[bool, str]:
        pass
