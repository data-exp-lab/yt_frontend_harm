"""Top-level package for yt frontend for HARM."""

__author__ = """Kacper Kowalik"""
__email__ = 'xarthisius.kk@gmail.com'
__version__ = '0.1.0'

from .data_structures import IHarmDataset, IHarmGrid, IHarmMKSHierarchy

from .fields import IHarmFieldInfo
from .io import IHarmIOHandler
from yt.utilities.object_registries import output_type_registry

add_harm_field = IHarmFieldInfo.add_field

output_type_registry["iharm"] = IHarmDataset
