OPENQASM 2.0;
include "qelib1.inc";
qreg q[50];
sdg q[18];
y q[29];
t q[39];
t q[15];
y q[5];
cx q[4],q[18];
swap q[18],q[45];
cz q[20],q[3];
cx q[47],q[37];
t q[49];
h q[3];
h q[13];
t q[28];
x q[13];
s q[19];
sdg q[17];
swap q[17],q[42];
swap q[36],q[38];
s q[2];
t q[34];
t q[18];
cx q[41],q[32];
t q[24];
x q[46];
cx q[4],q[28];
cz q[9],q[32];
swap q[43],q[17];
cx q[7],q[36];
swap q[2],q[16];
sdg q[12];
x q[36];
t q[20];
cx q[46],q[43];
z q[18];
swap q[11],q[33];
cz q[1],q[16];
t q[31];
cz q[43],q[11];
cz q[38],q[36];
cx q[20],q[19];
t q[45];
t q[39];
cz q[5],q[1];
t q[4];
cz q[37],q[35];
y q[4];
t q[9];
cz q[32],q[0];
y q[8];
cx q[22],q[31];
cz q[24],q[13];
swap q[8],q[7];
t q[48];
t q[19];
z q[9];
cz q[3],q[11];
t q[47];
x q[29];
s q[29];
x q[21];
z q[5];
t q[42];
t q[47];
t q[15];
t q[29];
cz q[28],q[31];
t q[28];
t q[38];
x q[5];
x q[37];
swap q[5],q[35];
z q[14];
t q[9];
t q[38];
sdg q[9];
swap q[9],q[35];
h q[43];
h q[45];
x q[37];
sdg q[12];
y q[40];
cz q[27],q[2];
sdg q[26];
cx q[48],q[7];
cz q[24],q[22];
z q[49];
t q[19];
t q[25];
swap q[13],q[19];
swap q[46],q[31];
sdg q[29];
cx q[42],q[15];
t q[45];
swap q[24],q[33];
h q[40];
h q[9];
t q[25];
z q[32];
cz q[35],q[46];
cz q[6],q[3];
y q[40];
t q[46];
t q[21];
cz q[35],q[12];
t q[34];
t q[1];
t q[17];
t q[25];
sdg q[2];
cx q[38],q[42];
h q[10];
h q[25];
x q[26];
t q[26];
t q[49];
t q[7];
t q[15];
t q[4];
x q[34];
t q[25];
z q[21];
h q[29];
h q[14];
t q[47];
z q[3];
z q[19];
x q[17];
t q[34];
y q[49];
t q[49];
t q[43];
swap q[34],q[20];
cz q[3],q[36];
swap q[11],q[26];
cz q[20],q[15];
cz q[28],q[32];
h q[25];
x q[15];
t q[24];
cx q[47],q[3];
t q[4];
h q[43];
cx q[29],q[35];
t q[33];
cx q[1],q[12];
t q[39];
swap q[29],q[34];
h q[49];
h q[41];
sdg q[45];
t q[0];
t q[9];
cx q[41],q[18];
t q[38];
swap q[41],q[25];
x q[40];
y q[30];
s q[49];
s q[27];
s q[29];
t q[1];
t q[42];
swap q[47],q[9];
t q[29];
h q[49];
swap q[6],q[4];
swap q[20],q[19];
sdg q[45];
swap q[11],q[16];
x q[11];
cz q[48],q[41];
t q[1];
cz q[12],q[46];
s q[34];
cz q[17],q[48];
t q[9];
t q[32];
s q[32];
sdg q[14];
cz q[49],q[46];
t q[29];
z q[7];
y q[30];
swap q[33],q[10];
swap q[27],q[10];
y q[31];
t q[15];
cz q[24],q[38];
cx q[36],q[19];
z q[27];
t q[26];
h q[25];
h q[22];
t q[42];
sdg q[32];
z q[28];
sdg q[14];
t q[8];
cz q[12],q[2];
t q[44];
y q[32];
y q[43];
t q[30];
z q[22];
s q[16];
x q[37];
h q[21];
s q[32];
h q[27];
cx q[0],q[4];
t q[22];
t q[37];
s q[20];
z q[25];
t q[14];
h q[10];
cz q[16],q[32];
h q[45];
t q[48];
t q[25];
cx q[9],q[16];
t q[46];
sdg q[20];
sdg q[34];
s q[23];
y q[6];
sdg q[14];
t q[31];
cx q[43],q[37];
swap q[46],q[38];
swap q[14],q[48];
cz q[27],q[49];
t q[40];
t q[30];
t q[30];
t q[26];
cz q[49],q[6];
z q[45];
cx q[36],q[33];
t q[49];
t q[48];
t q[27];
t q[4];
swap q[21],q[6];
swap q[12],q[4];
t q[34];
t q[24];
t q[13];
h q[6];
t q[43];
t q[10];
cx q[42],q[44];
cx q[49],q[16];
h q[45];
t q[46];
swap q[46],q[18];
cx q[8],q[41];
cx q[17],q[3];
t q[30];
t q[47];
z q[11];
cz q[47],q[42];
t q[31];
sdg q[41];
h q[6];
t q[22];
cx q[49],q[4];
t q[46];
swap q[3],q[28];
cz q[3],q[49];
t q[13];
y q[0];
sdg q[5];
s q[3];
cz q[41],q[30];
t q[19];
x q[26];
cz q[35],q[32];
z q[47];
t q[47];
s q[31];
swap q[27],q[13];
t q[22];
cz q[36],q[5];
y q[39];
t q[39];
t q[8];
sdg q[8];
t q[24];
t q[5];
t q[6];
cz q[4],q[14];
t q[25];
x q[34];
t q[48];
x q[6];
sdg q[20];
t q[23];
swap q[49],q[18];
t q[40];
t q[31];
cx q[2],q[3];
t q[9];
t q[18];
sdg q[6];
cz q[33],q[0];
t q[46];
t q[41];
cz q[1],q[38];
swap q[26],q[42];
cx q[7],q[6];
sdg q[7];
t q[42];
cz q[7],q[38];
h q[34];
cx q[1],q[5];
cz q[47],q[1];
z q[28];
swap q[6],q[13];
x q[3];
swap q[19],q[29];
t q[29];
s q[35];
z q[12];
t q[18];
x q[32];
t q[39];
h q[5];
h q[14];
t q[22];
cx q[26],q[48];
t q[42];
t q[1];
swap q[45],q[7];
t q[30];
swap q[23],q[40];
cx q[33],q[27];
s q[34];
x q[6];
t q[24];
t q[31];
swap q[49],q[48];
sdg q[18];
cz q[15],q[4];
t q[45];
swap q[6],q[47];
x q[31];
swap q[39],q[28];
swap q[36],q[4];
t q[49];
t q[46];
t q[0];
cz q[5],q[39];
sdg q[14];
cx q[22],q[26];
y q[8];
y q[14];
swap q[10],q[12];
x q[30];
y q[7];
sdg q[12];
cx q[27],q[13];
t q[37];
z q[17];
t q[0];
t q[32];
t q[15];
t q[30];
swap q[17],q[26];
swap q[36],q[44];
s q[31];
cx q[4],q[34];
swap q[2],q[12];
t q[24];
t q[3];
sdg q[21];
swap q[2],q[41];
swap q[33],q[45];
t q[47];
t q[26];
h q[38];
y q[29];
swap q[1],q[23];
cz q[24],q[11];
t q[24];
x q[28];
t q[26];
z q[26];
z q[31];
x q[21];
cx q[43],q[30];
x q[31];
t q[11];
x q[5];
cz q[20],q[28];
t q[18];
y q[4];
cz q[47],q[37];
swap q[36],q[6];
cx q[38],q[13];
t q[9];
cz q[25],q[6];
swap q[18],q[36];
t q[7];
x q[34];
t q[11];
t q[34];
t q[0];
t q[48];
cx q[36],q[10];
t q[27];
t q[37];
h q[33];
s q[13];
t q[45];
sdg q[17];
h q[35];
t q[5];
z q[47];
cz q[0],q[21];
s q[10];
t q[19];
cx q[19],q[27];
h q[14];
y q[34];
swap q[14],q[32];
z q[27];
t q[29];
cx q[0],q[46];
cx q[34],q[5];
cx q[44],q[15];
s q[48];
h q[43];
t q[37];
t q[26];
s q[30];
t q[17];
cz q[37],q[9];
swap q[4],q[9];
cz q[39],q[18];
cz q[34],q[44];
sdg q[38];
cx q[43],q[0];
t q[33];
s q[46];
sdg q[29];
s q[14];
swap q[42],q[41];
t q[31];
swap q[26],q[49];
h q[42];
cx q[49],q[35];
s q[19];
t q[46];
z q[7];
t q[21];
s q[45];
cx q[36],q[1];
swap q[13],q[6];
cx q[14],q[4];
cx q[5],q[4];
t q[36];
cz q[17],q[1];
x q[15];
x q[27];
swap q[38],q[17];
t q[9];
sdg q[1];
t q[39];
cz q[23],q[22];
cx q[31],q[37];
sdg q[5];
t q[33];
t q[33];
y q[28];
t q[43];
swap q[5],q[6];
cx q[34],q[17];
cz q[22],q[1];
cx q[4],q[49];
s q[18];
t q[5];
cz q[36],q[24];
t q[31];
swap q[19],q[44];
t q[28];
x q[45];
swap q[13],q[3];
t q[41];
y q[47];
s q[3];
t q[34];
z q[33];
cx q[30],q[15];
cx q[5],q[37];
cz q[48],q[36];
sdg q[38];
cx q[24],q[2];
cz q[33],q[36];
h q[11];
cx q[11],q[35];
h q[24];
t q[1];
t q[2];
y q[23];
cx q[39],q[42];
sdg q[14];
t q[28];
z q[41];
t q[14];
t q[37];
sdg q[16];
swap q[42],q[20];
cx q[24],q[17];
swap q[24],q[10];
t q[39];
t q[33];
cx q[18],q[46];
cx q[37],q[36];
cz q[49],q[41];
sdg q[41];
swap q[22],q[29];
sdg q[28];
h q[17];
h q[0];
swap q[15],q[4];
t q[33];
h q[25];
x q[45];
t q[4];
t q[35];
swap q[44],q[3];
t q[1];
t q[27];
y q[13];
cz q[1],q[37];
cz q[28],q[25];
t q[8];
h q[9];
swap q[26],q[6];
cx q[32],q[1];
h q[29];
y q[26];
sdg q[20];
t q[36];
cx q[42],q[5];
cx q[18],q[40];
h q[1];
t q[7];
h q[28];
z q[49];
swap q[20],q[9];
t q[8];
t q[6];
y q[22];
t q[2];
y q[23];
h q[44];
sdg q[33];
s q[5];
y q[30];
y q[3];
sdg q[29];
t q[15];
y q[45];
sdg q[28];
t q[34];
t q[9];
h q[27];
cz q[36],q[8];
h q[49];
swap q[37],q[40];
s q[0];
h q[4];
z q[19];
cz q[29],q[24];
cz q[49],q[10];
swap q[23],q[45];
cx q[26],q[7];
y q[14];
cx q[10],q[20];
cx q[15],q[21];
s q[35];
t q[45];
swap q[49],q[38];
s q[13];
t q[16];
z q[43];
t q[31];
y q[29];
sdg q[3];
s q[16];
cx q[35],q[19];
cx q[9],q[44];
swap q[33],q[13];
cx q[49],q[6];
x q[25];
t q[34];
s q[18];
y q[43];
t q[2];
x q[35];
h q[40];
t q[9];
sdg q[45];
y q[49];
cx q[24],q[14];
t q[5];
t q[16];
s q[40];
t q[19];
x q[45];
t q[49];
t q[34];
cx q[41],q[46];
swap q[37],q[1];
s q[27];
t q[49];
t q[18];
cz q[29],q[48];
cz q[22],q[44];
t q[5];
cz q[1],q[43];
t q[27];
swap q[4],q[46];
s q[1];
t q[2];
t q[45];
x q[1];
z q[48];
cx q[41],q[4];
swap q[40],q[15];
t q[31];
t q[39];
t q[19];
swap q[43],q[32];
y q[15];
swap q[32],q[31];
t q[30];
cx q[16],q[5];
y q[43];
cx q[30],q[5];
h q[1];
h q[2];
cx q[46],q[15];
cz q[36],q[29];
cz q[12],q[28];
y q[4];
cx q[24],q[17];
h q[18];
sdg q[9];
t q[9];
t q[8];
z q[6];
z q[23];
s q[8];
s q[40];
sdg q[45];
t q[44];
t q[14];
t q[32];
cz q[17],q[1];
sdg q[41];
y q[4];
z q[4];
cx q[38],q[21];
t q[35];
swap q[35],q[26];
sdg q[42];
cx q[38],q[41];
t q[27];
h q[27];
t q[44];
z q[21];
z q[41];
cz q[15],q[10];
t q[33];
s q[19];
cz q[23],q[13];
cz q[42],q[34];
h q[19];
cx q[32],q[0];
t q[13];
s q[49];
z q[21];
cz q[25],q[27];
swap q[39],q[12];
cz q[39],q[22];
y q[36];
z q[9];
y q[20];
t q[25];
t q[47];
z q[9];
t q[26];
t q[42];
t q[32];
cz q[5],q[36];
y q[45];
cx q[1],q[30];