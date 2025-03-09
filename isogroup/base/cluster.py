from typing import List, Union, Iterator, Self
from isogroup.base.feature import Feature
import re

class Cluster:

    def __init__(self, features: List[Feature], cluster_id = None, cluster_annotation = None, tracer = None, tracer_element = None):
        self.features = features
        self.cluster_id = cluster_id
        self.tracer_element = tracer_element
        self.tracer = tracer
        self.cluster_annotation = cluster_annotation # Name of the annotated cluster

    def __repr__(self) -> str:
        return f"Cluster({self.cluster_id}, {self.cluster_annotation}, {self.features})"
    
    
    def __len__(self) -> int:
        """
        Returns the number of unique features in the cluster, i.e the number of unique feature_id in the cluster
        """
        return len(set([f.feature_id for f in self.features]))


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
        Returns the metabolite objects associated to the annotated features in the cluster
        """
        return self.features[0].metabolite

    @property
    def f_annotations(self) -> str:
        """
        Returns the list of possible feature annotations in the cluster
        """
        return self.features[0].annotation
    
    @property
    def f_isotopologues(self) -> List[int]:
        """
        Returns the list of isotopologues for the annotated features in the cluster
        """
        return [f.isotopologue for f in self.features]

    @property
    def cluster_annotation_formula(self) -> str:
        """
        Returns the formula of the cluster annotation
        """
        # Check if the cluster is annotated
        if self.cluster_annotation is None:
            raise ValueError("No cluster annotation found. Run get_annotated_clusters() first")

        # Retrieve the formula for the metabolite matching the cluster_annotation
        for metabolite in self.metabolite:
            if metabolite.label == self.cluster_annotation:
                return metabolite.formula
        else:
            raise ValueError(f"No metabolite found with the label {self.cluster_annotation}")

    @property
    def cluster_annotation_isotopologues(self) -> List[int]:
        """
        Returns the list of isotopologues for feature.metabolite matching the cluster_annotation
        """
        cluster_isotopologues = set()

        for feature in self.features:
            if self.cluster_annotation in feature.annotation:
                for isotopologue_set in feature.isotopologue:
                    cluster_isotopologues.update(feature.isotopologue)
        return list(cluster_isotopologues)
                          

    @property
    def is_complete(self) -> bool:
        """
        Returns True if the cluster is complete
        """     
        formula = self.cluster_annotation_formula # Formula of the cluster annotation
        element_number = formula[self.tracer_element] # Number of tracer element in the formula
        
        return len(self) == element_number + 1

 
    @property
    def missing_isotopologue(self) -> List[int]:
        """
        Returns a list of missing isotopologues in the cluster
        """
        # Check if the cluster is complete
        if self.is_complete:
            return []
        
        formula = self.cluster_annotation_formula
        element_number = formula[self.tracer_element]
        expected_isotopologues = list(range(element_number + 1))
        current_isotopologues = self.cluster_annotation_isotopologues
        return [i for i in expected_isotopologues if i not in current_isotopologues]        
    
    @property
    def is_adduct(self) -> tuple[bool, str]:
        pass
