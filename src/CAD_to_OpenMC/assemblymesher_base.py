import abc

class assemblymesher(abc.ABC):

  verbosity_level=0

  @abc.abstractmethod
  def generate_stls():
    pass

  @classmethod
  def set_verbosity(cls,level:int=0):
    cls.verbosity_level=level
