from abc import ABC, abstractmethod


class Chain(ABC):
    @abstractmethod
    def __init__(self, protein_sequence: str) -> None:
        pass

    # `TODO
    # Discuss wheter this method is necessary in the Chain class
    @staticmethod
    @abstractmethod
    def build_turn_qubit() -> None:
        pass
