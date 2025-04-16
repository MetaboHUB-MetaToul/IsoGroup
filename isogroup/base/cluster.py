from typing import List, Union, Iterator, Self
from isogroup.base.feature import Feature

class Cluster:

    def __init__(self, features: list|None, cluster_id = None, name=None):
        self.features = features
        self.cluster_id = cluster_id
        self.tracer_element = features[0]._tracer_element if features is not None else None
        self.tracer = features[0].tracer if features is not None else None
        self.name = name
        self._formula = None

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

    def __remove__(self, other: Union[Feature, Self]) -> Self:
        """
        Remove a feature from the cluster
        :param other: Feature
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
    def metabolite(self):
        """
        Returns the complete annotation for each feature in the cluster
        """
        return [f.metabolite for f in self.features]
    
    @property
    def chemical(self):
        """
        Returns the chemical object of the cluster
        """
        for feature in self.features:
            if self.name in feature.metabolite:
                idx = feature.metabolite.index(self.name)
                return feature.chemical[idx]

    @property
    def element_number(self):
        """
        Returns the number of tracer element for the formula matching to the cluster name
        """
        self.formula
        return self._formula[self.tracer_element]

    
    @property
    def isotopologues(self) -> List[int]:
        """
        Returns the list of isotopologues in the cluster
        """
        isotopologues = []
        for feature in self.features:
            if self.name in feature.metabolite:
                idx = feature.metabolite.index(self.name)
                isotopologues.append(feature.isotopologue[idx])
        return isotopologues


    @property
    def formula(self) -> str:
        """
        Returns the formula of the metabolite matching to the cluster name
        """
        if self._formula is None:
            # Check if the cluster is annotated
            if self.name is None:
                raise ValueError("No cluster found. Run clusterize() first")
            
            if self.name in self.features[0].metabolite:
                idx = self.features[0].metabolite.index(self.name)
                self._formula = self.features[0].formula[idx]
            else:  
                raise ValueError(f"No metabolite matching '{self.name}' found in the cluster features")
            
        return self._formula


    @property
    def expected_isotopologues_in_cluster(self):
        """
        Returns the list of expected isotopologues in the cluster
        Based on the number of tracer element in its formula
        """
        return list(range(self.element_number + 1))
                           

    @property
    def is_complete(self) -> bool:
        """
        Returns True if the cluster is complete (i.e contains all isotopologues expected)
        """   
        return len(self) == self.element_number + 1 and self.expected_isotopologues_in_cluster == self.isotopologues
    
    @property
    def is_incomplete(self) -> bool:
        """
        Returns True if the cluster is incomplete (i.e contains not all isotopologues expected)
        """
        return len(self) < self.element_number + 1 or len(set(self.isotopologues)) != len(self.expected_isotopologues_in_cluster)

    @property
    def is_duplicated(self) -> bool:
        """
        Returns True if the cluster contains duplicated isotopologues
        """
        return len(set(self.isotopologues)) != len(self.isotopologues)

    @property
    def is_corrupted(self) -> bool:
        """
        Returns True if the cluster is corrupted (overfilled ?) (i.e contains isotopologues not expected)
        """ 
        pass


    @property 
    def status(self):
        """
        Returns the cluster status (Ok, incomplete, duplicated, corrupted)
        """
        if self.is_complete:
            return "Ok"
        
        if self.is_incomplete and self.is_duplicated:
            return "Incomplete, duplicated isotopologues"

        if self.is_incomplete:
            return "Incomplete"

        if self.is_duplicated:
            return "Duplicated isotopologues"
        
        # if self.is_corrupted:
        #     return "Corrupted"
        

    @property
    def missing_isotopologues(self) -> List[int]:
        """
        Returns a list of missing isotopologues in the cluster
        """
        if self.is_incomplete:
            return [i for i in self.expected_isotopologues_in_cluster if i not in self.isotopologues]
        

    @property
    def duplicated_isotopologues(self):
        """
        Returns a list of duplicated isotopologues in the cluster
        """
        if self.is_duplicated:
            return [i for i in set(self.isotopologues) if self.isotopologues.count(i) > 1]
        
    
    @property
    def is_adduct(self) -> tuple[bool, str]:
        pass

    @property
    def cluster_summary(self):
        """
        Returns a summary of the cluster
        """
        return {
            "cluster_id": self.cluster_id,
            "name": self.name,
            "number_of_features": len(self),
            "isotopologues": self.isotopologues,
            "status": self.status,
            "missing_isotopologues": self.missing_isotopologues,
            "duplicated_isotopologues": self.duplicated_isotopologues
        }