from qiskit import QuantumRegister, QuantumCircuit

chain = ["H", "P", "H", "P", "H"]
chain_length = len(chain)

configuration_qubits = 4 * (chain_length - 3)

configuration_register = QuantumRegister(configuration_qubits)

qc = QuantumCircuit(configuration_register)

for qubit in range(configuration_qubits):
    qc.h(qubit)

for qubit in range(configuration_qubits):
    qc.ry(0, qubit)

qc.barrier()

print(qc.draw())
