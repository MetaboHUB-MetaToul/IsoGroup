from typing import List, Union, Iterator, Self
from isogroup.base.feature import Feature

class Cluster:

    def __init__(self, features: list|None, cluster_id = None, name=None):
        self.features = features
        self.cluster_id = cluster_id
        self.tracer_element = features[0]._tracer_element if features is not None else None
        self.tracer = features[0].tracer if features is not None else None
        self.name = name
        self._cluster_formula = None

    def __repr__(self) -> str:
        return f"Cluster({self.cluster_id}, {self.features})"
    
    def __len__(self) -> int:
        """
        Returns the number of features in the cluster
        """
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

    @property
    def metabolite(self): ## NBU
        """
        Returns the complete annotation of each feature of the cluster
        """
        return [f.metabolite for f in self.features]
    
    @property
    def chemical(self):
        """
        Returns the chemical object of the cluster
        """
        return [f.chemical for f in self.features]

    @property
    def element_number(self):
        """
        Returns the number of tracer element in the cluster
        """
        self.cluster_formula
        return self._cluster_formula[self.tracer_element]

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
    def f_annotations(self) -> str:
        """
        Returns the list of possible feature annotations in the cluster
        """
        return self.features[0].metabolite
    
    @property
    def cluster_isotopologues(self) -> List[int]:
        """
        Returns the list of isotopologues for the annotated features in the cluster
        """
        isotopologues = []
        for feature in self.features:
            if self.name in feature.metabolite:
                idx = feature.metabolite.index(self.name)
                isotopologues.append(feature.isotopologue[idx])
        return isotopologues
        #return [f.isotopologue for f in self.features if f.metabolite == self.name]

    @property
    def cluster_formula(self) -> str:
        """
        Returns the feature.formula for the metabolite matching the cluster.name
        """
        if self._cluster_formula is None:
            # Check if the cluster is annotated
            if self.name is None:
                raise ValueError("No cluster found. Run clusterize() first")
            if len(self.features) == 0:
                raise ValueError("No features found in this cluster. Run clusterize() first")
            

            if self.name in self.features[0].metabolite:
                idx = self.features[0].metabolite.index(self.name)
                self._cluster_formula = self.features[0].formula[idx]
            else:  
                raise ValueError(f"No metabolite matching '{self.name}' found in the cluster features")
            
        return self._cluster_formula

    @property
    def expected_isotopologues_in_cluster(self):
        """
        Returns the list of expected isotopologues in the cluster
        """
        return list(range(self.element_number + 1))
                           

    @property
    def is_complete(self) -> bool:
        """
        Returns True if the cluster is complete (I.e contains all isotopologues expected)
        """   
        if len(self) == self.element_number + 1:
            if self.expected_isotopologues_in_cluster == self.cluster_isotopologues:
                return True
        return False
        #return len(self) == self.element_number + 1


    @property
    def missing_isotopologue(self) -> List[int]:
        """
        Returns a list of missing isotopologues in the cluster
        """
        # Check if the cluster is complete
        if self.is_complete:
            return []
        
        # expected_isotopologues = list(range(self.element_number + 1))
        # current_isotopologues = [f.isotopologue for f in self.features if f.metabolite == self.name]

        return [i for i in self.expected_isotopologues_in_cluster if i not in self.cluster_isotopologues]

    @property
    def duplicated_isotopologues(self):
        """
        Returns a list of duplicated isotopologues in the cluster
        """
        if self.is_complete:
            return []
        
        if self.expected_isotopologues_in_cluster == self.cluster_isotopologues:
            return []
        
        if self.expected_isotopologues_in_cluster != self.cluster_isotopologues:
            duplicated_isotopologues = [i for i in set(self.cluster_isotopologues) if self.cluster_isotopologues.count(i) > 1]
            # Return the list of duplicated isotopologues
            return duplicated_isotopologues
        
        # if self.is_complete:
        #     return []
        
        # expected_isotopologues = list(range(self.element_number + 1))
        # current_isotopologues = [f.isotopologue for f in self.features if f.metabolite == self.name]
        # return [i for i in expected_isotopologues if i not in current_isotopologues]        
    
    @property
    def is_adduct(self) -> tuple[bool, str]:
        pass


        # if self._cluster_metabolite_formula is None:
        #     # Check if the cluster is annotated
        #     if self.name is None:
        #         raise ValueError("No cluster metabolite found. Run clusterize() first")
        #     if len(self.features) == 0:
        #         raise ValueError("No features found in this cluster. Run clusterize() first")

        #     # Retrieve the formula for the metabolite matching the cluster_metabolite
        #     for feature in self.features:
        #         if feature.metabolite == self.name:
        #             return feature.formula
                
        # return self._cluster_metabolite_formula

            # @property
    # def cluster_annotation_isotopologues(self) -> List[int]:
    #     """
    #     Returns the list of isotopologues for feature.metabolite matching the cluster_annotation
    #     """
    #     cluster_isotopologues = set()

    #     for feature in self.features:
    #         if self.cluster_annotation in feature.annotation:
    #             for isotopologue_set in feature.isotopologue:
    #                 cluster_isotopologues.update(feature.isotopologue)
    #     return list(cluster_isotopologues)