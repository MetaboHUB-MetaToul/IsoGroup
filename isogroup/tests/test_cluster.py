from isogroup.base.cluster import Cluster
from isogroup.base.feature import Feature
import math

def test_lowest_highest_rt():
    """
    Test the lowest_rt and highest_rt properties of the Cluster class.
    """
    cluster = Cluster(features=[
        Feature(feature_id="F1", mz=119.02575, rt=3.5, intensity=1000, tracer="13C"),
        Feature(feature_id="F2", mz=120.02913, rt=1.2, intensity=500, tracer="13C"),
        Feature(feature_id="F3", mz=191.01958, rt=2.8, intensity=750, tracer="13C"),
    ], cluster_id="C1")
    assert cluster.lowest_rt == 1.2
    assert cluster.highest_rt == 3.5

def test_lowest_highest_mz():
    """
    Test the lowest_mz and highest_mz properties of the Cluster class.
    """
    cluster = Cluster(features=[
        Feature(feature_id="F1", mz=119.02575, rt=3.5, intensity=1000, tracer="13C"),
        Feature(feature_id="F2", mz=120.02913, rt=1.2, intensity=500, tracer="13C"),
        Feature(feature_id="F3", mz=191.01958, rt=2.8, intensity=750, tracer="13C"),
    ], cluster_id="C1")
    assert cluster.lowest_mz == 119.02575
    assert cluster.highest_mz == 191.01958

def test_mean_rt_mz():
    """
    Test the mean_rt and mean_mz properties of the Cluster class.
    """
    cluster = Cluster(features=[
        Feature(feature_id="F1", mz=119.02575, rt=3.5, intensity=1000, tracer="13C"),
        Feature(feature_id="F2", mz=120.02913, rt=1.2, intensity=500, tracer="13C"),
        Feature(feature_id="F3", mz=191.01958, rt=2.8, intensity=750, tracer="13C"),
    ], cluster_id="C1")

    assert math.isclose(cluster.mean_rt, (3.5 + 1.2 + 2.8) / 3)
    assert math.isclose(cluster.mean_mz, (119.02575 + 120.02913 + 191.01958) / 3)
