OPENQASM 2.0;
include "qelib1.inc";
qreg q[60];
t q[46];
cx q[35],q[57];
cz q[19],q[5];
swap q[31],q[26];
t q[21];
t q[30];
swap q[31],q[24];
h q[58];
t q[57];
y q[51];
cz q[25],q[15];
h q[3];
cz q[30],q[16];
t q[20];
t q[32];
cz q[9],q[28];
sdg q[26];
y q[15];
t q[51];
x q[50];
t q[15];
t q[47];
t q[43];
t q[56];
t q[33];
y q[50];
cz q[49],q[14];
s q[59];
x q[19];
cz q[30],q[25];
t q[49];
cx q[50],q[40];
t q[47];
t q[34];
swap q[3],q[53];
t q[3];
swap q[12],q[41];
cz q[36],q[5];
s q[19];
swap q[29],q[21];
s q[7];
t q[59];
h q[30];
h q[0];
cx q[52],q[50];
cx q[5],q[51];
sdg q[40];
sdg q[46];
cz q[11],q[39];
cx q[29],q[34];
t q[28];
sdg q[39];
t q[18];
cz q[25],q[4];
t q[14];
cx q[10],q[34];
sdg q[45];
h q[8];
t q[32];
t q[59];
h q[24];
x q[39];
t q[47];
s q[16];
x q[15];
cx q[45],q[44];
t q[55];
y q[19];
cz q[42],q[30];
h q[27];
t q[11];
swap q[37],q[31];
cx q[55],q[32];
t q[5];
swap q[5],q[46];
swap q[34],q[30];
z q[56];
t q[56];
y q[29];
y q[7];
swap q[12],q[14];
x q[31];
t q[7];
t q[10];
h q[37];
h q[46];
t q[26];
cz q[32],q[16];
sdg q[44];
cx q[13],q[7];
t q[12];
t q[0];
x q[54];
x q[51];
t q[24];
cx q[4],q[26];
cz q[37],q[39];
cz q[13],q[8];
cx q[37],q[55];
cx q[29],q[3];
swap q[2],q[54];
sdg q[17];
swap q[26],q[47];
swap q[53],q[23];
t q[52];
t q[40];
t q[0];
t q[54];
t q[3];
sdg q[17];
swap q[52],q[56];
sdg q[12];
s q[11];
h q[23];
t q[1];
z q[0];
t q[38];
swap q[54],q[28];
cx q[51],q[7];
t q[4];
swap q[33],q[31];
h q[56];
t q[57];
x q[28];
sdg q[34];
swap q[26],q[21];
t q[30];
h q[42];
x q[38];
swap q[13],q[18];
t q[48];
sdg q[22];
t q[3];
cx q[8],q[9];
z q[34];
s q[56];
h q[56];
x q[27];
swap q[13],q[59];
t q[32];
t q[35];
cx q[26],q[32];
cz q[34],q[11];
h q[28];
x q[24];
swap q[19],q[8];
cx q[43],q[1];
cx q[18],q[10];
cx q[40],q[5];
cz q[57],q[9];
t q[10];
t q[55];
t q[53];
cx q[58],q[43];
t q[36];
t q[34];
s q[15];
cz q[37],q[21];
t q[34];
y q[4];
swap q[41],q[50];
cx q[47],q[7];
cz q[14],q[7];
s q[20];
cx q[7],q[18];
t q[10];
t q[54];
cz q[46],q[23];
s q[10];
cz q[49],q[20];
s q[57];
h q[28];
z q[45];
t q[33];
t q[58];
t q[34];
h q[14];
swap q[46],q[23];
swap q[4],q[56];
t q[25];
t q[1];
t q[40];
cz q[42],q[24];
t q[47];
t q[40];
t q[7];
swap q[0],q[32];
y q[32];
t q[58];
t q[41];
cx q[22],q[3];
cx q[4],q[29];
t q[22];
swap q[57],q[46];
cx q[24],q[35];
t q[46];
s q[30];
sdg q[4];
t q[31];
t q[58];
cx q[52],q[15];
cz q[4],q[54];
cx q[19],q[5];
t q[41];
swap q[4],q[5];
swap q[26],q[28];
t q[45];
sdg q[25];
sdg q[44];
cx q[50],q[30];
t q[21];
cx q[54],q[36];
cx q[48],q[55];
t q[45];
swap q[4],q[17];
y q[26];
y q[51];
s q[53];
t q[10];
cx q[39],q[2];
cz q[24],q[22];
cz q[52],q[1];
cx q[49],q[59];
swap q[53],q[17];
cx q[14],q[39];
swap q[49],q[47];
swap q[49],q[26];
cx q[28],q[44];
t q[38];
t q[1];
swap q[3],q[58];
swap q[12],q[15];
t q[51];
t q[53];
y q[43];
x q[43];
cx q[7],q[58];
t q[17];
cx q[3],q[4];
z q[38];
s q[29];
x q[48];
x q[53];
t q[11];
t q[56];
cx q[41],q[18];
cx q[59],q[43];
cx q[43],q[12];
t q[56];
swap q[19],q[22];
t q[22];
swap q[54],q[53];
swap q[12],q[45];
z q[27];
x q[3];
cz q[56],q[27];
cx q[19],q[53];
t q[39];
t q[25];
x q[22];
t q[6];
x q[43];
t q[32];
s q[55];
swap q[3],q[40];
h q[3];
swap q[50],q[25];
sdg q[26];
swap q[9],q[27];
t q[30];
swap q[27],q[49];
swap q[21],q[48];
cx q[49],q[44];
t q[30];
s q[49];
t q[24];
y q[51];
t q[33];
t q[15];
swap q[35],q[31];
sdg q[1];
t q[20];
t q[19];
z q[8];
cz q[31],q[38];
sdg q[52];
cx q[58],q[25];
swap q[2],q[25];
y q[30];
y q[50];
swap q[52],q[55];
y q[10];
h q[25];
y q[43];
cx q[13],q[55];
swap q[13],q[54];
t q[13];
t q[51];
y q[50];
s q[2];
z q[37];
s q[48];
cx q[57],q[35];
t q[14];
sdg q[29];
cx q[38],q[58];
sdg q[21];
t q[25];
cz q[33],q[45];
t q[16];
t q[3];
t q[34];
z q[27];
sdg q[52];
t q[42];
t q[43];
x q[28];
t q[0];
t q[25];
s q[57];
sdg q[8];
x q[7];
t q[4];
s q[28];
s q[11];
swap q[30],q[36];
t q[33];
cz q[14],q[45];
sdg q[1];
y q[19];
t q[53];
t q[44];
swap q[49],q[36];
cz q[5],q[22];
h q[0];
s q[32];
t q[23];
h q[44];
swap q[15],q[40];
swap q[33],q[52];
cz q[12],q[48];
t q[38];
h q[47];
y q[25];
t q[43];
h q[25];
x q[57];
t q[21];
cx q[50],q[19];
t q[48];
t q[19];
cx q[58],q[3];
cx q[38],q[48];
cx q[30],q[50];
t q[42];
y q[3];
swap q[53],q[26];
t q[58];
cz q[13],q[6];
cz q[11],q[50];
s q[36];
z q[26];
sdg q[28];
t q[33];
cz q[52],q[36];
t q[8];
cx q[26],q[29];
swap q[8],q[10];
swap q[42],q[50];
t q[59];
t q[1];
s q[15];
t q[55];
t q[1];
t q[58];
z q[23];
s q[5];
x q[1];
t q[56];
t q[0];
cz q[16],q[51];
z q[15];
x q[10];
cx q[38],q[5];
cz q[17],q[44];
swap q[46],q[38];
t q[6];
swap q[33],q[21];
t q[14];
sdg q[54];
t q[49];
h q[23];
cz q[6],q[47];
t q[51];
cz q[12],q[44];
h q[31];
cz q[38],q[0];
h q[37];
t q[59];
t q[59];
swap q[32],q[25];
z q[3];
t q[48];
cz q[28],q[31];
swap q[15],q[43];
t q[29];
swap q[23],q[17];
sdg q[26];
t q[12];
z q[56];
cx q[20],q[16];
t q[17];
t q[2];
t q[51];
t q[5];
swap q[50],q[0];
cz q[42],q[50];
cx q[28],q[36];
cz q[5],q[44];
sdg q[9];
t q[17];
cz q[7],q[21];
swap q[7],q[11];
cz q[18],q[16];
cz q[55],q[12];
t q[57];
t q[31];
h q[35];
swap q[10],q[39];
t q[43];
t q[29];
cx q[17],q[25];
z q[35];
swap q[31],q[37];
s q[10];
cx q[41],q[55];
cx q[7],q[4];
t q[31];
swap q[23],q[55];
t q[52];
t q[19];
cx q[31],q[56];
t q[3];
h q[50];
t q[8];
x q[3];
y q[45];
t q[50];
t q[35];
cz q[2],q[9];
z q[44];
h q[0];
h q[20];
cx q[7],q[12];
t q[53];
sdg q[50];
h q[24];
swap q[49],q[14];
swap q[4],q[34];
t q[52];
s q[13];
cz q[46],q[24];
cx q[31],q[11];
t q[14];
t q[58];
sdg q[3];
h q[42];
t q[40];
cx q[30],q[3];
cz q[51],q[59];
t q[55];
t q[1];
sdg q[31];
t q[37];
cz q[36],q[29];
cz q[2],q[35];
t q[1];
cx q[34],q[42];
t q[41];
t q[7];
swap q[15],q[6];
cx q[27],q[50];
t q[3];
s q[48];
s q[2];
swap q[7],q[2];
cx q[52],q[38];
t q[36];
x q[48];
h q[9];
t q[54];
t q[44];
t q[36];
z q[23];
t q[12];
z q[15];
sdg q[18];
cz q[3],q[52];
swap q[8],q[48];
t q[18];
t q[25];
cx q[44],q[9];
cx q[9],q[29];
t q[50];
h q[10];
cz q[29],q[4];
swap q[31],q[46];
cz q[49],q[26];
cz q[31],q[48];
cx q[48],q[30];
y q[50];
cz q[13],q[40];
s q[30];
sdg q[47];
t q[17];
swap q[39],q[54];
swap q[9],q[28];
s q[13];
sdg q[38];
t q[24];
t q[24];
swap q[29],q[36];
h q[8];
swap q[52],q[36];
swap q[0],q[25];
x q[30];
t q[2];
t q[34];
cx q[3],q[8];
t q[23];
cx q[5],q[47];
swap q[23],q[46];
cx q[51],q[29];
z q[15];
cz q[26],q[18];
t q[27];
cz q[55],q[23];
t q[6];
s q[26];
cz q[5],q[37];
h q[47];
h q[20];
cx q[23],q[47];
swap q[25],q[21];
t q[37];
cx q[41],q[0];
y q[9];
y q[56];
y q[46];
y q[45];
cx q[9],q[24];
s q[50];
t q[53];
sdg q[23];
y q[11];
t q[48];
swap q[31],q[46];
cz q[57],q[39];
t q[54];
t q[25];
t q[19];
cx q[51],q[40];
t q[51];
t q[55];
t q[59];
h q[52];
sdg q[8];
z q[11];
swap q[32],q[5];
swap q[0],q[1];
cz q[53],q[18];
t q[40];
z q[20];
z q[53];
x q[56];
t q[27];
swap q[20],q[12];
swap q[11],q[42];
cx q[21],q[58];
cz q[17],q[7];
s q[32];
cx q[32],q[50];
cx q[28],q[58];
s q[26];
t q[33];
swap q[28],q[17];
cx q[34],q[56];
t q[48];
cz q[43],q[13];
sdg q[4];
swap q[42],q[45];
t q[6];
sdg q[15];
cz q[41],q[27];
t q[11];
t q[34];
cx q[48],q[25];
cx q[44],q[25];
cz q[0],q[57];
sdg q[26];
cz q[49],q[58];
swap q[28],q[33];
swap q[11],q[18];
h q[36];
cz q[43],q[30];
t q[53];
sdg q[58];
t q[40];
y q[26];
swap q[8],q[40];
h q[7];
t q[32];
cx q[7],q[10];
z q[21];
t q[59];
t q[52];
swap q[36],q[27];
x q[34];
cx q[35],q[5];
t q[2];
t q[28];
t q[32];
swap q[54],q[24];
swap q[49],q[30];
swap q[11],q[52];
cz q[46],q[20];
t q[36];
t q[12];
cz q[0],q[49];
y q[30];
cz q[48],q[21];
sdg q[31];
swap q[11],q[14];
t q[39];
swap q[26],q[44];
x q[33];
h q[15];
s q[16];
cz q[23],q[38];
cz q[36],q[44];
cz q[50],q[37];
z q[3];
y q[26];
swap q[50],q[0];
t q[21];
t q[1];
x q[12];
t q[47];
z q[3];
z q[3];
cz q[33],q[3];
cx q[7],q[24];
t q[14];
cz q[19],q[55];
h q[32];
s q[9];
s q[49];
cz q[12],q[34];
cx q[21],q[7];
s q[43];
t q[51];
x q[41];
t q[38];
t q[37];
swap q[40],q[30];
swap q[22],q[40];
t q[42];
cx q[43],q[45];
cz q[46],q[8];
sdg q[27];
cx q[21],q[14];
t q[31];
swap q[45],q[37];
swap q[37],q[33];
t q[11];
cx q[40],q[23];
swap q[32],q[43];
t q[45];
swap q[45],q[52];
x q[15];
swap q[4],q[49];
s q[6];
t q[48];
t q[28];
y q[16];
h q[57];
t q[19];
y q[21];
swap q[49],q[15];
h q[42];
sdg q[4];
s q[43];
cz q[11],q[24];
z q[46];
s q[54];
t q[24];
t q[0];
t q[6];
swap q[44],q[45];
cx q[37],q[25];
t q[1];
h q[41];
t q[33];
cx q[59],q[58];
cx q[44],q[30];
t q[3];
x q[27];
cx q[40],q[21];
cx q[6],q[30];
cz q[22],q[18];
t q[45];
cx q[22],q[51];
h q[35];
t q[2];
t q[51];
y q[42];
t q[18];
t q[40];
t q[25];
z q[9];
h q[46];
sdg q[46];
sdg q[26];
cx q[7],q[51];
x q[16];
t q[35];
s q[40];
t q[44];
t q[25];
t q[14];
t q[31];
s q[3];
y q[53];
sdg q[47];
t q[19];
cz q[28],q[48];
t q[37];
t q[17];
t q[2];
sdg q[5];
cx q[40],q[37];
s q[7];
cz q[50],q[37];
t q[46];
swap q[3],q[35];
z q[12];
swap q[40],q[10];
cz q[52],q[58];
x q[47];
t q[35];
z q[59];
h q[55];
cz q[19],q[5];
t q[53];
x q[53];
sdg q[27];
y q[13];
x q[52];
y q[14];
t q[2];
t q[52];
t q[9];
cz q[51],q[31];
t q[30];
z q[15];
t q[28];
sdg q[1];
t q[51];
t q[18];
cx q[4],q[46];
cz q[30],q[31];
t q[19];
t q[42];
cz q[11],q[57];
swap q[22],q[58];
cz q[33],q[19];
cx q[22],q[53];
h q[16];
z q[55];
swap q[26],q[47];
t q[57];
t q[28];
h q[40];
s q[41];
t q[3];
t q[1];
x q[33];
y q[48];
t q[25];
t q[50];
t q[0];
cx q[48],q[26];
t q[48];
t q[26];
t q[53];
y q[0];
x q[52];
cz q[12],q[0];
z q[15];
t q[17];
t q[40];
x q[59];
s q[25];
t q[13];
h q[0];
cx q[13],q[49];
cx q[31],q[4];
x q[33];
t q[28];
h q[0];
t q[48];
t q[4];
s q[12];
y q[48];
z q[20];
sdg q[10];
cz q[18],q[43];
cx q[53],q[1];
swap q[26],q[0];
y q[47];
swap q[22],q[30];
cz q[1],q[2];
swap q[23],q[43];
cx q[35],q[7];
t q[34];
t q[12];
swap q[1],q[25];
swap q[3],q[40];
sdg q[16];
z q[20];
cx q[42],q[35];
cz q[39],q[8];
t q[8];
cx q[27],q[40];
h q[13];
t q[42];
cx q[40],q[7];
s q[53];
t q[1];
cx q[25],q[56];
cz q[18],q[47];
swap q[33],q[44];
swap q[17],q[55];
t q[20];
t q[15];
h q[55];
t q[22];
swap q[20],q[27];
cz q[1],q[17];
t q[16];
swap q[39],q[54];
t q[48];
cx q[57],q[53];
cz q[19],q[16];
t q[33];
y q[14];
sdg q[5];
x q[11];
t q[29];
t q[39];
cx q[0],q[22];
h q[4];
cx q[29],q[32];
t q[15];
t q[42];
t q[59];
z q[8];
t q[43];
t q[47];
swap q[4],q[55];
swap q[7],q[56];
t q[30];
cx q[43],q[34];
s q[19];
z q[49];
x q[6];
y q[40];
x q[15];
t q[9];
swap q[26],q[42];
cz q[18],q[58];
t q[45];
z q[25];
y q[12];
cx q[55],q[11];
t q[40];
t q[39];
t q[36];
cz q[27],q[51];
s q[40];
swap q[32],q[10];
x q[2];
s q[27];
t q[18];
swap q[47],q[3];
sdg q[37];
h q[51];
h q[14];
cz q[51],q[48];
z q[31];