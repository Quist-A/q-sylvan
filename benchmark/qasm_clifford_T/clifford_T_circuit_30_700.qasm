OPENQASM 2.0;
include "qelib1.inc";
qreg q[30];
t q[19];
swap q[18],q[2];
cx q[17],q[27];
cz q[19],q[16];
y q[0];
swap q[5],q[13];
t q[16];
y q[29];
swap q[16],q[29];
t q[17];
swap q[0],q[23];
t q[14];
x q[18];
swap q[14],q[17];
z q[9];
sdg q[24];
h q[4];
x q[28];
t q[8];
cz q[16],q[22];
swap q[22],q[19];
s q[26];
s q[15];
t q[27];
x q[26];
h q[4];
cz q[7],q[11];
cz q[29],q[15];
z q[19];
swap q[29],q[5];
t q[10];
cx q[6],q[1];
t q[8];
t q[7];
x q[25];
t q[24];
cz q[15],q[8];
t q[2];
z q[16];
cx q[29],q[23];
x q[8];
x q[9];
t q[20];
t q[28];
swap q[24],q[12];
s q[8];
t q[6];
cx q[0],q[28];
t q[9];
cz q[29],q[2];
sdg q[3];
x q[17];
cx q[24],q[29];
t q[27];
t q[13];
t q[20];
t q[23];
t q[9];
swap q[13],q[7];
t q[24];
y q[21];
swap q[17],q[27];
s q[23];
s q[23];
t q[20];
swap q[21],q[7];
cz q[12],q[2];
s q[29];
cx q[4],q[13];
t q[4];
t q[2];
t q[18];
x q[25];
z q[8];
cz q[7],q[24];
t q[13];
t q[29];
t q[16];
swap q[13],q[27];
t q[17];
swap q[16],q[19];
cz q[19],q[0];
s q[7];
t q[18];
h q[16];
swap q[12],q[1];
t q[22];
cx q[15],q[23];
swap q[24],q[0];
t q[10];
t q[9];
s q[23];
cx q[3],q[21];
y q[19];
x q[9];
t q[6];
t q[22];
t q[3];
cx q[27],q[26];
t q[3];
y q[10];
t q[5];
cx q[24],q[28];
t q[22];
sdg q[19];
cz q[7],q[6];
x q[18];
t q[18];
cz q[17],q[27];
t q[14];
h q[24];
cz q[6],q[4];
swap q[19],q[3];
sdg q[1];
x q[0];
t q[20];
s q[10];
cz q[2],q[8];
z q[27];
cx q[0],q[17];
y q[19];
y q[24];
t q[28];
h q[12];
cx q[15],q[1];
cx q[4],q[14];
s q[6];
t q[21];
h q[28];
x q[18];
swap q[3],q[11];
cz q[12],q[1];
z q[11];
h q[8];
cx q[27],q[5];
cx q[13],q[24];
swap q[22],q[15];
t q[22];
cx q[12],q[14];
s q[8];
sdg q[27];
h q[20];
x q[20];
t q[26];
z q[16];
swap q[27],q[20];
cz q[12],q[10];
s q[26];
y q[5];
cx q[3],q[27];
x q[29];
sdg q[20];
cz q[25],q[7];
h q[11];
cz q[27],q[5];
swap q[11],q[10];
x q[5];
cx q[7],q[4];
cx q[9],q[13];
s q[28];
sdg q[20];
z q[22];
h q[29];
t q[17];
sdg q[28];
cz q[10],q[29];
cx q[24],q[29];
y q[22];
cz q[2],q[15];
x q[15];
cx q[28],q[27];
s q[26];
z q[29];
t q[23];
t q[25];
cx q[25],q[10];
cz q[4],q[16];
h q[11];
t q[27];
y q[8];
t q[12];
swap q[2],q[9];
y q[17];
swap q[26],q[25];
t q[18];
cx q[18],q[21];
cx q[21],q[13];
cx q[2],q[12];
cz q[18],q[22];
z q[7];
t q[23];
cx q[10],q[3];
swap q[28],q[23];
z q[19];
t q[28];
sdg q[17];
t q[4];
swap q[4],q[11];
t q[27];
z q[13];
sdg q[16];
t q[6];
s q[10];
y q[0];
swap q[17],q[2];
t q[24];
cx q[7],q[3];
swap q[6],q[9];
t q[29];
t q[6];
t q[15];
cx q[7],q[16];
h q[11];
t q[18];
t q[5];
cz q[18],q[12];
cz q[16],q[21];
t q[22];
t q[19];
x q[7];
x q[8];
cz q[14],q[17];
t q[21];
t q[26];
h q[21];
swap q[29],q[12];
cx q[27],q[14];
cx q[16],q[18];
h q[11];
cx q[27],q[0];
t q[0];
sdg q[18];
t q[10];
x q[19];
z q[24];
z q[22];
s q[13];
t q[0];
x q[18];
t q[26];
t q[18];
cx q[6],q[8];
cx q[27],q[9];
t q[16];
t q[10];
t q[27];
t q[7];
h q[22];
swap q[5],q[9];
h q[3];
t q[10];
t q[11];
swap q[13],q[15];
x q[12];
s q[15];
t q[5];
swap q[4],q[14];
t q[21];
swap q[12],q[11];
sdg q[5];
cx q[22],q[12];
x q[7];
t q[23];
h q[23];
t q[15];
z q[22];
cz q[23],q[6];
t q[1];
z q[11];
cx q[28],q[10];
t q[20];
t q[2];
t q[25];
cx q[16],q[15];
cx q[6],q[11];
y q[27];
swap q[8],q[11];
swap q[21],q[1];
cx q[26],q[12];
sdg q[13];
z q[3];
sdg q[4];
s q[21];
cx q[12],q[2];
t q[26];
t q[6];
cx q[22],q[19];
t q[14];
cx q[12],q[13];
t q[18];
h q[29];
x q[9];
t q[20];
cz q[8],q[17];
swap q[3],q[6];
t q[1];
t q[13];
swap q[0],q[26];
t q[11];
t q[2];
t q[15];
swap q[23],q[25];
t q[20];
s q[7];
y q[14];
t q[19];
cz q[2],q[27];
h q[14];
y q[24];
t q[12];
t q[11];
t q[19];
t q[25];
t q[7];
x q[14];
t q[5];
t q[19];
cz q[23],q[2];
swap q[4],q[25];
cx q[5],q[29];
t q[17];
h q[0];
x q[12];
cz q[9],q[23];
t q[24];
cx q[15],q[17];
swap q[6],q[2];
t q[11];
h q[17];
cx q[6],q[9];
h q[24];
cx q[27],q[6];
cz q[5],q[20];
t q[11];
x q[16];
cz q[3],q[26];
s q[11];
cx q[0],q[2];
t q[9];
swap q[17],q[16];
t q[15];
t q[27];
t q[22];
t q[10];
cx q[13],q[9];
cz q[5],q[8];
t q[15];
y q[22];
y q[20];
s q[21];
cz q[2],q[19];
t q[10];
t q[18];
cx q[24],q[9];
t q[29];
swap q[22],q[6];
cx q[24],q[28];
h q[10];
cz q[20],q[17];
t q[11];
t q[16];
cz q[21],q[6];
swap q[28],q[1];
sdg q[8];
cz q[27],q[0];
y q[25];
cz q[20],q[3];
s q[4];
h q[23];
cz q[18],q[5];
cz q[26],q[4];
z q[4];
t q[6];
s q[6];
cz q[9],q[23];
swap q[5],q[15];
t q[29];
t q[16];
t q[17];
swap q[9],q[1];
t q[4];
h q[9];
x q[0];
t q[17];
h q[13];
t q[12];
t q[17];
t q[28];
sdg q[2];
s q[14];
t q[5];
sdg q[0];
t q[14];
x q[4];
cz q[19],q[26];
h q[4];
x q[16];
cz q[5],q[24];
t q[9];
y q[17];
t q[7];
t q[29];
swap q[20],q[17];
t q[3];
t q[10];
cx q[5],q[10];
t q[15];
s q[15];
cz q[11],q[4];
cx q[5],q[24];
x q[29];
y q[14];
y q[17];
t q[4];
s q[1];
swap q[26],q[14];
sdg q[10];
cx q[2],q[27];
t q[17];
swap q[22],q[5];
s q[29];
swap q[19],q[13];
t q[15];
swap q[23],q[26];
t q[4];
t q[14];
s q[1];
cz q[17],q[11];
t q[17];
cx q[24],q[8];
cx q[4],q[26];
swap q[9],q[10];
x q[3];
t q[26];
z q[24];
swap q[14],q[16];
h q[0];
sdg q[0];
z q[8];
t q[12];
cz q[14],q[20];
t q[14];
h q[7];
cx q[12],q[2];
cz q[4],q[12];
x q[8];
cz q[22],q[17];
cz q[5],q[4];
cx q[23],q[28];
swap q[7],q[22];
swap q[24],q[27];
swap q[22],q[3];
cz q[9],q[12];
swap q[26],q[25];
s q[14];
t q[3];
cz q[16],q[13];
sdg q[12];
x q[16];
cz q[10],q[7];
s q[8];
cx q[6],q[7];
t q[7];
swap q[8],q[16];
cx q[12],q[4];
cx q[12],q[5];
h q[19];
y q[12];
swap q[2],q[6];
y q[9];
t q[29];
t q[29];
cx q[29],q[24];
s q[8];
sdg q[19];
t q[11];
t q[9];
z q[1];
swap q[23],q[20];
z q[14];
h q[12];
cx q[5],q[10];
z q[23];
t q[16];
t q[0];
cx q[2],q[15];
t q[9];
cz q[16],q[0];
h q[25];
cx q[23],q[13];
cz q[10],q[8];
cx q[27],q[28];
swap q[13],q[26];
cx q[29],q[0];
cz q[26],q[14];
cz q[21],q[0];
t q[19];
swap q[16],q[29];
cx q[21],q[2];
t q[2];
cx q[5],q[17];
x q[22];
t q[21];
t q[11];
z q[18];
cx q[13],q[1];
t q[11];
t q[13];
t q[21];
cx q[28],q[27];
swap q[23],q[18];
cx q[27],q[20];
sdg q[14];
t q[29];
t q[18];
t q[29];
cx q[3],q[17];
swap q[18],q[27];
cx q[10],q[2];
cz q[25],q[11];
t q[22];
t q[25];
t q[22];
y q[27];
swap q[3],q[22];
s q[22];
sdg q[3];
swap q[23],q[29];
t q[15];
y q[13];
s q[16];
y q[16];
sdg q[6];
h q[29];
t q[27];
sdg q[19];
t q[13];
cx q[14],q[26];
sdg q[21];
t q[1];
h q[25];
cx q[14],q[2];
cx q[7],q[26];
t q[23];
x q[1];
t q[18];
cz q[13],q[16];
swap q[3],q[24];
y q[8];
cz q[29],q[3];
z q[15];
sdg q[10];
cx q[13],q[0];
cx q[29],q[5];
s q[11];
s q[21];
t q[6];
y q[3];
y q[12];
h q[6];
x q[9];
h q[6];
swap q[1],q[2];
t q[10];
cx q[27],q[8];
cz q[7],q[1];
cz q[5],q[11];
z q[28];
sdg q[11];
swap q[2],q[20];
t q[13];
x q[20];
t q[1];
sdg q[18];
cz q[9],q[20];
cx q[20],q[23];
t q[15];
t q[0];
t q[29];
s q[29];
z q[27];
t q[13];
t q[15];
h q[14];
t q[10];
cz q[11],q[16];
x q[10];
x q[5];
swap q[28],q[9];
sdg q[20];
t q[19];
t q[22];
cx q[10],q[6];
t q[0];
x q[26];
cx q[5],q[11];
cz q[13],q[17];
sdg q[7];
y q[11];
t q[20];
sdg q[15];
swap q[7],q[10];
t q[22];
cz q[26],q[22];
z q[18];
z q[1];
t q[20];
cx q[20],q[21];
t q[2];
t q[22];
x q[28];
cz q[17],q[16];
x q[9];
t q[28];
sdg q[23];
sdg q[29];
y q[1];
sdg q[8];
h q[0];
sdg q[2];
t q[7];
t q[16];
y q[14];
t q[24];
y q[11];
cz q[27],q[6];
s q[21];
swap q[0],q[6];
cz q[5],q[1];
t q[19];
cz q[20],q[21];
z q[20];
x q[1];
cz q[17],q[27];
h q[17];
t q[7];
t q[16];
cz q[16],q[22];
cz q[29],q[0];
z q[7];
cx q[14],q[18];
t q[18];
h q[20];
s q[11];
cx q[23],q[0];
y q[7];
cx q[23],q[0];
s q[23];
t q[21];
cx q[20],q[7];
t q[17];
t q[12];
z q[29];
h q[12];
t q[25];
s q[15];
cx q[25],q[0];
cz q[4],q[10];
t q[23];
swap q[23],q[10];
cz q[9],q[5];
cx q[17],q[14];
y q[23];
y q[27];
swap q[17],q[6];
sdg q[18];
h q[20];
h q[11];
t q[20];
s q[14];
t q[19];
cx q[3],q[1];
t q[27];
s q[29];
t q[14];
t q[14];
s q[7];
t q[7];
z q[10];
cx q[12],q[14];
x q[0];
y q[18];
swap q[8],q[4];
t q[12];
x q[7];
t q[26];
t q[22];
sdg q[12];
swap q[28],q[6];
t q[6];
cz q[19],q[16];
cz q[14],q[5];
t q[10];
cx q[25],q[17];
h q[17];
t q[12];
cx q[20],q[25];
cx q[22],q[21];
swap q[11],q[20];
t q[14];