OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
cz q[2],q[14];
cz q[19],q[6];
x q[9];
cz q[4],q[18];
cz q[1],q[14];
x q[5];
t q[10];
y q[15];
z q[1];
x q[3];
z q[17];
t q[11];
t q[1];
swap q[10],q[12];
cz q[11],q[3];
t q[16];
t q[16];
cx q[18],q[14];
sdg q[14];
t q[3];
t q[4];
x q[19];
t q[10];
cz q[0],q[3];
x q[12];
cx q[6],q[9];
swap q[15],q[14];
t q[4];
cx q[7],q[1];
cx q[4],q[17];
s q[9];
y q[19];
cz q[11],q[7];
x q[6];
t q[4];
t q[16];
z q[0];
z q[5];
cz q[15],q[4];
cx q[2],q[12];
t q[18];
t q[1];
t q[14];
s q[19];
t q[6];
x q[7];
cx q[12],q[5];
t q[14];
cz q[3],q[14];
t q[11];
t q[12];
sdg q[16];
cz q[7],q[1];
sdg q[11];
t q[4];
cx q[0],q[16];
y q[4];
t q[19];
t q[19];
h q[17];
cx q[14],q[6];
cz q[5],q[2];
y q[2];
cz q[6],q[3];
y q[10];
sdg q[17];
swap q[10],q[18];
t q[3];
cz q[13],q[15];
z q[16];
t q[12];
cz q[1],q[7];
x q[19];
t q[17];
t q[2];
h q[6];
cz q[13],q[8];
s q[12];
cx q[6],q[11];
cz q[1],q[10];
s q[10];
sdg q[14];
z q[4];
h q[16];
t q[9];
swap q[9],q[11];
cz q[8],q[9];
cx q[11],q[7];
sdg q[7];
h q[6];
swap q[0],q[16];
t q[2];
z q[4];
cz q[17],q[4];
s q[2];
x q[15];
x q[3];
t q[10];
t q[3];
t q[1];
y q[13];
cz q[9],q[0];
t q[5];
s q[0];
sdg q[18];
z q[2];
sdg q[0];
swap q[0],q[13];
s q[6];
y q[1];
t q[16];
t q[3];
cx q[8],q[17];
t q[19];
x q[6];
t q[9];
cz q[4],q[8];
t q[16];
cx q[1],q[8];
x q[6];
t q[10];
t q[12];
t q[5];
x q[8];
t q[1];
t q[8];
swap q[17],q[18];
h q[0];
y q[14];
cz q[12],q[3];
z q[17];
x q[8];
cz q[2],q[0];
s q[18];
swap q[7],q[17];
t q[15];
z q[3];
h q[1];
z q[4];
swap q[16],q[5];
t q[17];
t q[19];
cz q[2],q[11];
t q[6];
swap q[9],q[1];
t q[6];
h q[0];
s q[18];
cz q[5],q[17];
cx q[1],q[11];
sdg q[16];
h q[10];
t q[18];
t q[4];
swap q[2],q[3];
s q[0];
h q[16];
cx q[14],q[3];
t q[9];
cx q[8],q[6];
cx q[8],q[15];
swap q[2],q[3];
t q[12];
y q[13];
cz q[14],q[13];
swap q[2],q[1];
cz q[0],q[15];
s q[2];
t q[18];
t q[8];
cz q[9],q[15];
cz q[19],q[3];
cz q[0],q[2];
t q[17];
cx q[8],q[10];
t q[16];
cx q[1],q[6];
x q[15];
y q[17];
cz q[3],q[19];
swap q[3],q[1];
cz q[18],q[17];
cx q[19],q[14];
t q[7];
t q[2];
s q[2];
cz q[11],q[0];
cx q[1],q[8];
cz q[9],q[15];
swap q[19],q[13];
x q[4];
swap q[0],q[8];
t q[8];
t q[10];
x q[13];
swap q[19],q[6];
z q[4];
t q[12];
y q[12];
swap q[5],q[12];
y q[9];
cz q[15],q[10];
s q[11];
t q[0];
cz q[1],q[0];
t q[4];
t q[7];
cz q[17],q[3];
cz q[8],q[4];
t q[2];
cx q[6],q[11];
x q[10];
t q[17];
x q[18];
swap q[15],q[5];
t q[0];
swap q[3],q[12];
y q[15];
t q[15];
swap q[5],q[14];
t q[15];
h q[12];
s q[2];
t q[15];
swap q[10],q[13];
cx q[7],q[9];
t q[19];
y q[12];
swap q[9],q[14];
h q[5];
y q[2];
x q[5];
h q[15];
t q[9];
t q[10];
cx q[8],q[10];
y q[12];
y q[19];
s q[6];
t q[19];
t q[14];
cx q[17],q[9];
cx q[4],q[10];
t q[1];
t q[18];
cx q[17],q[6];
cz q[18],q[14];
s q[15];
h q[13];
cz q[9],q[17];
t q[0];
t q[1];
s q[5];
swap q[19],q[5];
t q[2];
y q[11];
t q[13];
y q[16];
t q[11];
sdg q[2];
t q[18];
swap q[9],q[17];
t q[1];
h q[12];
t q[14];
t q[16];
t q[8];
t q[13];
t q[9];
sdg q[15];
t q[14];
cz q[17],q[15];
sdg q[10];
t q[3];
swap q[4],q[19];
cx q[18],q[14];
x q[8];
swap q[17],q[5];
t q[7];
cz q[0],q[5];
cz q[7],q[12];
t q[3];
swap q[11],q[17];
t q[4];
cz q[6],q[19];
t q[6];
cz q[5],q[13];
t q[12];
cx q[13],q[11];
sdg q[5];
cz q[17],q[1];
cz q[10],q[18];
x q[9];
s q[2];
t q[12];
cz q[16],q[4];
swap q[4],q[13];
swap q[12],q[1];
cx q[12],q[10];
cx q[16],q[2];
cz q[6],q[19];
swap q[10],q[2];
t q[19];
y q[10];
t q[13];
sdg q[18];
cz q[10],q[12];
t q[15];
s q[17];
z q[3];
z q[15];
cx q[10],q[7];
y q[16];
x q[10];
cx q[0],q[11];
t q[12];
sdg q[1];
cx q[5],q[0];
sdg q[7];
h q[9];
cx q[13],q[16];
swap q[7],q[12];
t q[10];
t q[3];
cz q[10],q[7];
h q[13];
t q[1];
t q[7];
t q[18];
h q[8];
t q[7];
x q[16];
t q[12];
s q[8];
z q[14];
s q[9];
cx q[6],q[7];
t q[14];
cz q[12],q[8];
s q[11];
x q[17];
s q[18];
cx q[1],q[7];
swap q[3],q[15];
cz q[13],q[6];
t q[11];
cx q[18],q[0];
s q[5];
y q[19];
sdg q[11];
cz q[14],q[8];
y q[14];
s q[8];
t q[3];
t q[11];
t q[0];
cz q[0],q[15];
x q[5];
t q[19];
h q[6];
cz q[15],q[8];
cz q[2],q[0];
swap q[14],q[3];
swap q[6],q[14];
t q[10];
swap q[8],q[9];
cz q[3],q[14];
cx q[17],q[11];
swap q[1],q[10];
swap q[0],q[11];
s q[1];
t q[3];
swap q[2],q[19];
sdg q[0];
y q[17];
swap q[1],q[14];
cx q[4],q[9];
swap q[11],q[6];
t q[18];
t q[14];
sdg q[3];
h q[0];
z q[8];
t q[14];
sdg q[13];
t q[9];
s q[11];
t q[8];
cz q[4],q[8];
t q[10];
t q[11];
t q[2];
h q[14];
cz q[16],q[13];
y q[9];
swap q[6],q[3];
t q[5];
t q[14];
y q[3];
t q[5];
h q[19];
sdg q[11];
t q[3];
sdg q[4];
t q[12];
x q[1];
s q[16];
swap q[7],q[14];
swap q[15],q[17];
t q[8];
swap q[4],q[16];
t q[8];
h q[6];
x q[3];
z q[15];
cz q[5],q[15];
x q[8];
cz q[16],q[3];
cx q[15],q[3];
y q[17];
cz q[7],q[2];
cx q[7],q[4];
h q[9];
t q[12];
cx q[18],q[12];
h q[10];
t q[9];
swap q[1],q[0];
t q[12];
t q[15];
t q[16];
t q[15];
z q[6];
t q[11];
t q[10];
cz q[14],q[13];
sdg q[9];
swap q[10],q[14];
s q[13];
t q[15];
cz q[3],q[4];
sdg q[11];
s q[10];
h q[15];
h q[4];
x q[10];
sdg q[2];
t q[8];
t q[0];
s q[3];
cz q[19],q[9];
t q[4];
s q[11];
y q[17];
h q[13];
swap q[16],q[12];
t q[10];
z q[16];
swap q[5],q[15];
t q[11];
cx q[16],q[12];
t q[0];
t q[1];
z q[9];
cx q[13],q[4];
cz q[15],q[3];
cz q[6],q[14];
t q[19];
cx q[0],q[2];
swap q[6],q[16];
cz q[1],q[18];
swap q[18],q[7];
t q[9];
t q[13];
y q[10];
swap q[5],q[18];
cz q[16],q[3];
t q[3];
t q[17];
t q[16];
sdg q[5];
t q[14];
swap q[19],q[14];
sdg q[13];
sdg q[17];
h q[18];
cz q[0],q[7];
swap q[1],q[0];
s q[14];
t q[19];
t q[12];
y q[2];
t q[15];
t q[11];
cz q[17],q[4];
t q[3];
t q[8];
sdg q[0];
t q[5];
swap q[0],q[3];
cx q[13],q[2];
t q[9];
t q[14];
swap q[8],q[10];
s q[7];
cx q[19],q[17];
cz q[15],q[0];
cx q[16],q[13];
t q[9];
t q[19];
cx q[8],q[9];
cx q[5],q[0];
swap q[5],q[0];
y q[7];
x q[13];
x q[9];
cz q[1],q[8];
t q[13];
s q[7];
cz q[1],q[15];
z q[17];
cz q[12],q[6];
y q[3];
y q[5];
t q[18];
t q[5];
t q[15];
cz q[8],q[12];
t q[14];
s q[18];
t q[6];
t q[1];
sdg q[13];
cx q[14],q[2];
t q[3];
cx q[8],q[17];
sdg q[6];
x q[7];
t q[6];
t q[4];
h q[0];
s q[11];
swap q[4],q[11];
cx q[5],q[19];
swap q[18],q[12];
h q[3];
swap q[4],q[6];
t q[10];
z q[19];
t q[11];
x q[17];
cz q[0],q[12];
t q[10];
h q[7];
t q[13];
cz q[4],q[8];
t q[9];
s q[5];
cz q[14],q[0];
cx q[18],q[13];
z q[11];
t q[9];
sdg q[8];
t q[11];
cx q[17],q[18];
t q[13];
t q[11];
t q[8];
t q[14];
h q[15];
t q[19];
sdg q[12];
s q[0];
cz q[19],q[6];
cz q[17],q[2];
swap q[13],q[2];
t q[9];
z q[1];
t q[18];
y q[5];
t q[4];
z q[4];
cz q[5],q[4];
x q[2];
h q[2];
cx q[0],q[16];
t q[17];
t q[17];
cx q[14],q[15];
z q[14];
t q[1];
x q[3];
y q[10];
x q[14];
t q[1];
cz q[17],q[3];
sdg q[4];
t q[18];
t q[14];
t q[7];
sdg q[19];
swap q[7],q[18];
h q[11];
h q[13];
cz q[1],q[0];
h q[13];
t q[9];
t q[16];
swap q[12],q[3];
h q[6];
t q[19];
cz q[0],q[2];
cx q[2],q[6];
cx q[10],q[16];
t q[16];
z q[9];
s q[17];
y q[2];
t q[1];
cx q[0],q[1];
t q[9];
cz q[6],q[8];
x q[18];
swap q[14],q[8];
h q[2];
t q[16];
z q[17];
swap q[6],q[18];
t q[12];
x q[5];
swap q[3],q[8];
cx q[9],q[6];
cx q[0],q[4];
y q[16];
t q[10];
cx q[18],q[17];
swap q[14],q[19];
t q[4];
cx q[10],q[16];
t q[4];
z q[3];
t q[11];
x q[5];
t q[18];
t q[9];
t q[0];
swap q[3],q[4];
cz q[5],q[16];
cz q[12],q[10];
t q[19];
cz q[17],q[3];
t q[14];
swap q[11],q[19];
sdg q[3];
cx q[19],q[6];
z q[18];
cz q[17],q[10];
cx q[14],q[5];
cz q[13],q[14];
swap q[1],q[9];
t q[16];
y q[9];
t q[0];
t q[18];
t q[16];
swap q[16],q[12];
y q[15];
swap q[2],q[10];
x q[5];
x q[11];
h q[3];
t q[17];
s q[13];
swap q[15],q[5];
t q[3];
swap q[15],q[0];
t q[10];
swap q[7],q[4];
t q[1];
x q[0];
h q[18];
h q[11];
cz q[10],q[18];
y q[18];
t q[10];
swap q[10],q[13];
t q[11];
y q[8];
t q[8];
swap q[12],q[18];
t q[16];
t q[0];
swap q[19],q[0];
t q[3];
cz q[16],q[1];
t q[7];
sdg q[14];
sdg q[15];
t q[5];
t q[18];
cz q[6],q[8];
t q[10];
h q[4];
z q[13];
swap q[14],q[16];
x q[1];
t q[4];
cz q[19],q[10];
cz q[17],q[0];
cx q[4],q[0];
cz q[10],q[16];
t q[4];
t q[12];
s q[19];
s q[6];
y q[2];
swap q[10],q[7];
t q[9];
s q[12];
t q[3];
sdg q[19];
swap q[17],q[1];
h q[12];
cz q[1],q[13];
h q[17];
t q[0];
swap q[5],q[15];
sdg q[7];
t q[11];
t q[15];
x q[5];
h q[11];
h q[12];
s q[13];
t q[3];
swap q[17],q[14];
t q[15];
y q[17];
sdg q[0];
cz q[14],q[10];
y q[8];
cx q[8],q[11];
t q[1];
z q[18];
swap q[19],q[12];
t q[15];
sdg q[7];
t q[1];
t q[17];
y q[9];
t q[9];
h q[6];
s q[0];
cx q[5],q[17];
swap q[0],q[15];
swap q[8],q[18];
s q[8];
t q[11];
t q[16];
y q[9];
s q[3];
swap q[5],q[13];
t q[15];
s q[16];
cz q[19],q[16];
y q[9];
t q[3];
t q[14];
z q[18];
x q[19];
t q[4];
x q[16];
swap q[5],q[17];
swap q[8],q[2];
h q[12];
y q[17];
cz q[3],q[9];
cz q[1],q[16];
x q[6];
sdg q[9];
cx q[19],q[1];
cx q[0],q[7];
t q[16];
cx q[0],q[18];
y q[12];
swap q[13],q[9];
t q[2];
t q[5];
swap q[16],q[5];
swap q[19],q[17];
t q[19];
cx q[4],q[13];
cx q[15],q[14];
sdg q[9];
sdg q[6];
sdg q[0];
cz q[16],q[1];
cx q[1],q[5];
t q[2];
z q[1];
swap q[17],q[6];
t q[12];
t q[7];
z q[5];
t q[9];
y q[15];
h q[2];
z q[4];
z q[18];
cx q[15],q[12];
z q[3];
cz q[3],q[7];
z q[11];
y q[16];
swap q[7],q[19];
x q[6];
z q[11];
s q[7];
z q[5];
t q[14];
s q[19];
s q[6];
y q[1];
z q[13];
sdg q[9];
y q[13];
t q[4];
swap q[10],q[4];
cz q[5],q[1];
t q[14];
t q[16];
cx q[13],q[18];
s q[6];
t q[9];
cz q[4],q[12];
h q[1];
t q[2];
y q[19];
cx q[5],q[3];
sdg q[3];
cz q[18],q[9];
t q[10];
t q[16];
cx q[11],q[18];
t q[0];
cx q[2],q[18];
cx q[0],q[13];
z q[15];
z q[15];
t q[5];
y q[7];
x q[8];
cz q[4],q[2];
sdg q[17];
x q[4];
cz q[16],q[15];
cz q[5],q[2];
t q[16];
sdg q[6];
swap q[5],q[18];
t q[13];
z q[16];
swap q[10],q[0];
cx q[18],q[1];
cz q[12],q[10];
sdg q[4];
cx q[2],q[9];
cx q[0],q[18];
t q[2];
t q[14];
t q[11];
t q[19];
sdg q[18];
swap q[19],q[3];
cz q[15],q[6];
t q[16];
s q[5];
y q[14];
cx q[4],q[9];
t q[12];
t q[15];
x q[17];
t q[8];
cx q[13],q[0];
s q[3];
s q[16];
swap q[15],q[12];
cx q[15],q[9];
y q[9];
z q[12];
cx q[7],q[14];
h q[13];
cz q[4],q[2];
cz q[18],q[10];
cx q[16],q[10];
t q[17];
cx q[8],q[17];
cx q[19],q[1];
t q[11];
cz q[17],q[0];
h q[8];
t q[5];
t q[1];
t q[9];
sdg q[19];
t q[11];
s q[11];
sdg q[13];
h q[8];
h q[6];
h q[2];
cz q[11],q[0];
s q[3];
x q[1];
h q[5];
t q[16];
x q[5];
cz q[11],q[2];
s q[10];
sdg q[16];
t q[19];
y q[8];
t q[19];
t q[1];
y q[9];
t q[17];
sdg q[11];
t q[10];
cx q[8],q[13];
t q[18];
cz q[1],q[19];
s q[1];
cx q[14],q[18];
t q[12];
cx q[14],q[6];
h q[5];
y q[16];
t q[6];
cx q[17],q[2];
cz q[0],q[5];
swap q[8],q[16];
t q[0];
t q[17];
t q[6];
y q[17];
cz q[14],q[0];
t q[0];
t q[17];
t q[6];
t q[19];