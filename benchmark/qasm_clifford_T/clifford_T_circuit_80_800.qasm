OPENQASM 2.0;
include "qelib1.inc";
qreg q[80];
t q[62];
t q[48];
cz q[41],q[74];
y q[45];
x q[33];
cx q[74],q[41];
z q[23];
swap q[7],q[55];
cx q[74],q[77];
s q[79];
swap q[30],q[4];
t q[21];
z q[46];
t q[64];
cx q[64],q[41];
t q[67];
t q[28];
y q[31];
t q[5];
swap q[23],q[36];
swap q[40],q[61];
cz q[2],q[34];
cz q[21],q[14];
t q[63];
cx q[19],q[75];
cz q[18],q[40];
cz q[62],q[59];
swap q[60],q[6];
t q[32];
swap q[30],q[4];
cz q[52],q[79];
cz q[49],q[7];
y q[65];
t q[28];
t q[55];
swap q[15],q[68];
swap q[31],q[58];
t q[75];
s q[31];
t q[42];
cx q[26],q[59];
h q[42];
cz q[54],q[72];
cz q[23],q[9];
cx q[67],q[48];
sdg q[12];
cz q[57],q[42];
sdg q[12];
x q[49];
cz q[16],q[46];
cx q[11],q[26];
t q[19];
x q[27];
cx q[37],q[40];
x q[23];
y q[62];
cz q[63],q[75];
s q[14];
t q[45];
t q[2];
cx q[4],q[26];
cx q[12],q[23];
swap q[44],q[7];
t q[2];
h q[45];
t q[67];
cz q[11],q[62];
cx q[26],q[3];
x q[8];
t q[9];
s q[13];
z q[33];
y q[31];
cz q[24],q[17];
z q[20];
t q[31];
sdg q[45];
sdg q[11];
y q[6];
swap q[75],q[50];
cx q[29],q[48];
swap q[60],q[33];
t q[9];
s q[52];
s q[60];
z q[29];
swap q[29],q[70];
t q[38];
y q[38];
t q[48];
t q[68];
cz q[50],q[31];
x q[12];
t q[0];
t q[10];
cx q[27],q[63];
h q[40];
z q[42];
t q[57];
cx q[45],q[56];
cx q[77],q[53];
h q[65];
swap q[27],q[6];
z q[56];
x q[73];
cz q[41],q[27];
t q[72];
swap q[21],q[18];
t q[79];
swap q[21],q[72];
cx q[47],q[16];
t q[54];
t q[18];
cz q[25],q[77];
y q[2];
swap q[42],q[74];
swap q[54],q[44];
h q[23];
t q[46];
cz q[15],q[10];
cz q[33],q[40];
cz q[61],q[38];
swap q[36],q[60];
cz q[62],q[0];
s q[32];
y q[15];
cz q[74],q[52];
sdg q[48];
t q[18];
swap q[31],q[22];
y q[57];
s q[18];
t q[50];
h q[50];
z q[1];
swap q[13],q[54];
cz q[43],q[61];
y q[59];
t q[59];
cz q[40],q[78];
sdg q[58];
y q[45];
t q[8];
swap q[24],q[49];
y q[6];
swap q[10],q[7];
t q[19];
sdg q[61];
sdg q[33];
cx q[65],q[66];
cx q[76],q[55];
swap q[10],q[58];
t q[31];
x q[35];
t q[41];
cz q[1],q[68];
cz q[18],q[60];
t q[70];
cz q[4],q[49];
t q[66];
sdg q[45];
t q[31];
t q[22];
y q[22];
cx q[75],q[10];
cz q[26],q[60];
t q[53];
t q[50];
z q[16];
cz q[2],q[52];
cx q[32],q[24];
cx q[4],q[3];
t q[30];
swap q[25],q[71];
swap q[34],q[38];
t q[79];
swap q[32],q[61];
sdg q[69];
cz q[25],q[51];
t q[24];
cz q[20],q[61];
t q[56];
swap q[23],q[71];
swap q[12],q[17];
swap q[40],q[29];
sdg q[2];
t q[63];
swap q[57],q[55];
t q[55];
h q[14];
y q[34];
t q[23];
cz q[11],q[55];
t q[9];
x q[64];
cz q[37],q[53];
y q[39];
cx q[36],q[6];
cx q[66],q[18];
swap q[25],q[10];
cz q[72],q[74];
cx q[36],q[76];
t q[64];
cz q[24],q[69];
t q[14];
s q[73];
t q[45];
t q[36];
swap q[16],q[77];
cz q[25],q[79];
swap q[25],q[24];
z q[37];
t q[16];
cz q[17],q[15];
cx q[38],q[55];
z q[79];
swap q[58],q[5];
t q[3];
t q[2];
swap q[20],q[43];
t q[37];
s q[14];
cx q[17],q[28];
t q[36];
h q[29];
swap q[67],q[9];
s q[58];
t q[1];
swap q[53],q[57];
s q[24];
sdg q[24];
h q[68];
swap q[43],q[14];
swap q[36],q[29];
cx q[19],q[18];
swap q[60],q[24];
swap q[5],q[11];
t q[38];
swap q[68],q[50];
cz q[16],q[6];
y q[14];
swap q[49],q[57];
cz q[64],q[71];
t q[46];
t q[78];
h q[66];
h q[58];
t q[0];
swap q[25],q[54];
y q[39];
cx q[14],q[48];
t q[40];
t q[44];
cz q[51],q[52];
t q[39];
z q[23];
cx q[33],q[24];
cx q[41],q[59];
t q[69];
z q[62];
cx q[32],q[69];
s q[24];
t q[73];
cx q[9],q[35];
z q[56];
cz q[76],q[26];
cx q[25],q[63];
swap q[1],q[31];
cz q[40],q[61];
s q[52];
t q[25];
t q[79];
s q[63];
sdg q[3];
sdg q[50];
swap q[76],q[66];
t q[3];
y q[5];
sdg q[63];
t q[7];
s q[68];
t q[8];
cx q[75],q[77];
swap q[18],q[43];
t q[10];
cz q[32],q[57];
t q[46];
s q[34];
t q[39];
s q[32];
cx q[45],q[79];
t q[68];
t q[43];
h q[60];
cz q[8],q[58];
swap q[10],q[21];
s q[14];
cx q[29],q[60];
t q[79];
sdg q[24];
t q[63];
t q[50];
s q[76];
t q[52];
t q[51];
t q[56];
t q[2];
cx q[10],q[39];
s q[70];
s q[16];
x q[32];
z q[78];
t q[4];
swap q[46],q[76];
t q[25];
y q[22];
t q[47];
y q[18];
y q[53];
h q[31];
t q[46];
cx q[72],q[67];
t q[69];
t q[5];
cx q[66],q[20];
swap q[12],q[6];
t q[41];
t q[73];
swap q[29],q[64];
t q[71];
cz q[27],q[26];
cz q[46],q[65];
cx q[30],q[18];
swap q[41],q[32];
t q[47];
cz q[7],q[23];
cz q[46],q[31];
h q[31];
s q[35];
z q[73];
t q[13];
t q[32];
t q[22];
s q[41];
x q[53];
z q[79];
cz q[33],q[54];
swap q[47],q[44];
cx q[69],q[51];
x q[70];
swap q[67],q[16];
t q[4];
t q[63];
x q[42];
y q[20];
cx q[0],q[18];
t q[26];
cz q[24],q[4];
cx q[18],q[72];
t q[65];
s q[42];
t q[3];
s q[8];
cz q[74],q[11];
t q[0];
t q[64];
swap q[31],q[56];
sdg q[53];
t q[62];
cx q[42],q[64];
h q[10];
cz q[15],q[3];
swap q[27],q[53];
x q[41];
t q[28];
cx q[31],q[38];
t q[78];
swap q[55],q[19];
z q[78];
cz q[56],q[4];
t q[42];
h q[66];
t q[74];
s q[18];
cz q[41],q[15];
y q[22];
s q[64];
t q[45];
cz q[1],q[27];
t q[67];
t q[43];
t q[39];
t q[16];
t q[5];
cx q[49],q[48];
t q[79];
cz q[47],q[58];
cx q[74],q[10];
t q[4];
t q[53];
cx q[61],q[71];
t q[77];
t q[64];
cz q[22],q[23];
cz q[59],q[70];
x q[47];
s q[48];
cx q[79],q[47];
t q[36];
t q[7];
cz q[1],q[34];
x q[41];
cz q[38],q[75];
t q[10];
t q[74];
cx q[61],q[72];
swap q[62],q[4];
t q[26];
t q[10];
t q[62];
t q[35];
y q[63];
t q[37];
swap q[22],q[75];
z q[6];
cz q[64],q[7];
s q[44];
x q[17];
sdg q[31];
swap q[66],q[65];
sdg q[40];
y q[70];
t q[28];
t q[75];
t q[71];
t q[53];
t q[23];
sdg q[49];
swap q[15],q[39];
t q[56];
cx q[49],q[68];
swap q[66],q[36];
t q[78];
t q[7];
x q[21];
cx q[53],q[79];
z q[9];
cx q[23],q[60];
cx q[33],q[67];
cz q[19],q[37];
t q[53];
t q[5];
s q[66];
y q[60];
swap q[13],q[75];
s q[56];
t q[7];
cz q[32],q[44];
cz q[34],q[63];
cx q[70],q[31];
h q[15];
cx q[36],q[61];
cx q[2],q[12];
t q[24];
t q[25];
swap q[41],q[62];
z q[63];
sdg q[17];
sdg q[49];
swap q[5],q[61];
swap q[18],q[64];
x q[55];
h q[2];
cz q[42],q[16];
h q[75];
cx q[5],q[0];
h q[23];
h q[73];
t q[75];
t q[12];
swap q[63],q[32];
t q[25];
cz q[55],q[69];
z q[2];
cx q[38],q[10];
swap q[70],q[49];
cx q[65],q[71];
sdg q[66];
cx q[38],q[9];
cx q[44],q[53];
swap q[3],q[59];
t q[51];
t q[27];
cx q[7],q[19];
swap q[61],q[10];
swap q[13],q[69];
z q[45];
swap q[12],q[30];
x q[66];
y q[49];
swap q[75],q[13];
s q[5];
t q[47];
t q[74];
t q[66];
cz q[52],q[43];
z q[54];
s q[31];
z q[72];
swap q[19],q[79];
h q[59];
swap q[47],q[16];
cz q[51],q[35];
swap q[57],q[35];
t q[51];
y q[39];
swap q[54],q[79];
s q[9];
t q[6];
t q[67];
t q[74];
z q[64];
sdg q[64];
cx q[62],q[51];
y q[63];
t q[34];
t q[46];
z q[37];
h q[70];
swap q[28],q[15];
x q[34];
cx q[77],q[42];
swap q[0],q[33];
cx q[36],q[54];
swap q[13],q[37];
t q[66];
sdg q[51];
t q[55];
cx q[44],q[14];
t q[6];
cx q[43],q[75];
t q[36];
cx q[9],q[22];
t q[2];
swap q[30],q[25];
t q[16];
cx q[17],q[76];
t q[45];
cz q[52],q[78];
cz q[68],q[6];
s q[67];
swap q[44],q[46];
h q[49];
sdg q[69];
t q[53];
cx q[27],q[51];
cz q[72],q[73];
t q[7];
cz q[6],q[36];
swap q[29],q[6];
z q[55];
cx q[20],q[64];
swap q[7],q[44];
cz q[75],q[69];
t q[63];
t q[56];
swap q[63],q[73];
t q[58];
cx q[46],q[18];
cz q[16],q[77];
cz q[14],q[18];
cx q[59],q[16];
t q[30];
h q[6];
s q[16];
s q[70];
s q[34];
s q[9];
t q[78];
z q[28];
s q[57];
t q[40];
cx q[32],q[64];
t q[8];
s q[6];
swap q[21],q[74];
t q[30];
sdg q[10];
t q[44];
t q[77];
cz q[54],q[40];
cx q[31],q[5];
t q[70];
t q[33];
cz q[15],q[23];
cx q[37],q[77];
cz q[60],q[32];
t q[29];
y q[31];
t q[28];
cz q[7],q[8];
t q[48];
z q[17];
y q[27];
z q[18];
cz q[37],q[26];
z q[17];
swap q[5],q[59];
cz q[76],q[58];
s q[79];
t q[32];
z q[43];
z q[42];
swap q[2],q[21];
y q[61];
cx q[36],q[35];
z q[41];
t q[32];
y q[67];
swap q[36],q[13];
swap q[15],q[56];
t q[42];
t q[71];
t q[5];
t q[7];
sdg q[52];
y q[49];
t q[46];
cz q[75],q[54];
h q[57];
x q[59];
t q[72];
t q[62];
h q[4];
cz q[6],q[22];
t q[13];
cz q[48],q[78];
x q[69];
t q[0];
cx q[21],q[13];
x q[21];
t q[58];
cx q[37],q[78];
y q[62];
h q[68];
t q[26];
cz q[45],q[65];
t q[50];
t q[73];
h q[9];
t q[30];
t q[21];
t q[25];
t q[37];
t q[6];
s q[77];
t q[68];
cx q[2],q[30];
sdg q[10];
swap q[43],q[58];
z q[51];
x q[30];
t q[13];
x q[74];
t q[12];
t q[53];
t q[36];
z q[34];
x q[67];
cx q[15],q[59];
x q[54];
t q[6];
t q[77];
swap q[9],q[65];
s q[61];
swap q[69],q[40];
t q[0];
swap q[70],q[9];
t q[69];
t q[38];
z q[13];
x q[71];
t q[22];
y q[67];
swap q[76],q[4];
t q[59];
t q[30];
y q[50];
y q[72];
cx q[51],q[56];
x q[21];
y q[26];
t q[75];
t q[37];
z q[17];
t q[63];
s q[51];
t q[77];
swap q[7],q[57];
t q[49];
t q[15];
sdg q[32];
y q[24];
t q[8];
s q[13];
cz q[25],q[27];
t q[16];
cz q[38],q[49];
cx q[1],q[70];
y q[31];
cz q[38],q[76];
sdg q[47];
cx q[52],q[73];
sdg q[33];
s q[25];
h q[53];
t q[61];
swap q[28],q[32];
y q[44];
t q[51];
t q[74];
cx q[25],q[59];
y q[47];
x q[20];
cz q[24],q[79];
t q[28];
y q[18];
cz q[27],q[39];
cz q[47],q[65];
cz q[33],q[62];
x q[0];
t q[43];
s q[74];
z q[54];
t q[44];
sdg q[9];
t q[47];
z q[29];
t q[25];
t q[45];
h q[76];
cz q[41],q[36];
s q[25];
t q[78];
t q[47];
y q[18];
t q[23];
t q[65];
t q[36];
sdg q[51];
cx q[54],q[43];
s q[76];
z q[66];
z q[24];
swap q[26],q[64];
x q[39];
z q[42];
x q[77];
cx q[67],q[42];
cx q[25],q[34];
z q[65];
y q[57];
t q[25];
t q[61];
swap q[6],q[14];
t q[27];
t q[16];
swap q[14],q[8];
cz q[62],q[25];
s q[44];
y q[69];
t q[41];
cx q[54],q[29];
t q[74];
cz q[73],q[10];
swap q[24],q[19];
cz q[20],q[79];
sdg q[24];
t q[66];
t q[62];
t q[13];
cx q[51],q[61];
t q[56];
y q[45];
t q[15];
x q[54];
t q[40];
t q[32];
cx q[37],q[7];
cz q[44],q[41];
y q[5];
t q[71];
t q[74];
z q[33];
t q[72];
x q[75];
cx q[67],q[31];
t q[61];
t q[2];
t q[34];