from abc import ABC, abstractmethod, abstractproperty

class Experiment(ABC):
    def __init__(self) -> None:
        self._conditions = {}
    
    @property
    def conditions(self):
        return self._conditions
    
    @conditions.setter
    def conditions(self, condition_data: dict):
        self._conditions = condition_data

    @abstractmethod
    def run_condition(self, condition, condition_data):
        ...

    def run(self, callback = None):
        result = {}
        for condition, condition_data in self._conditions.items():
            result[condition] = self.run_condition(condition, condition_data)
            if callback is not None:
                callback(result[condition])

        return result
