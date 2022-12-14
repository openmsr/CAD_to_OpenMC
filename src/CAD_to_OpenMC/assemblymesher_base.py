import abc

class assemblymesher(abc.ABC):
  @abc.abstractmethod
  def generate_stls():
    pass

  @classmethod
  def set_verbosity(cls,level):
    cls.verbosity_level=level
