#!/bin/bash
for j in {14..20}
do
    num_gates=$((j*50))
    for i in {2..8}
    do
        num_qubits=$((i*10))
        echo "Simulating num_qubits = " ${num_qubits} "num_gates = " ${num_gates}
        timeout 60 ../build/qasm/run_qasm_on_qmdd qasm_clifford_T/clifford_T_circuit_${num_qubits}_${num_gates}.qasm -j benchmark_results/${num_qubits}_${num_gates}_float_low-omega.json -c -e float -s low
        echo "Finished  low     float "
        timeout 60 ../build/qasm/run_qasm_on_qmdd qasm_clifford_T/clifford_T_circuit_${num_qubits}_${num_gates}.qasm -j benchmark_results/${num_qubits}_${num_gates}_float_l2-omega.json -c -e float -s l2
        echo "Finished  l2      float "
        timeout 60 ../build/qasm/run_qasm_on_qmdd qasm_clifford_T/clifford_T_circuit_${num_qubits}_${num_gates}.qasm -j benchmark_results/${num_qubits}_${num_gates}_qisq2_low-omega.json -c -e qisq2 -s low
        echo "Finished  low     qisq2 "
        echo "Finished simulation num_qubits = " ${num_qubits} "num_gates = " ${num_gates}
    done
done
