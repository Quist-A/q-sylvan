OPENQASM 2.0;
include "qelib1.inc";
qreg q[40];
x q[16];
cx q[3],q[39];
t q[25];
x q[28];
t q[26];
s q[28];
y q[0];
z q[19];
z q[38];
cz q[1],q[20];
t q[15];
cz q[34],q[31];
s q[6];
cz q[27],q[6];
x q[16];
cz q[31],q[2];
cz q[33],q[15];
x q[24];
t q[10];
swap q[24],q[11];
t q[39];
h q[15];
y q[11];
t q[7];
cz q[37],q[7];
z q[35];
t q[20];
swap q[4],q[27];
cz q[3],q[7];
swap q[39],q[33];
t q[4];
s q[16];
h q[17];
h q[1];
h q[18];
t q[15];
t q[4];
t q[6];
t q[11];
cz q[30],q[25];
y q[0];
x q[13];
cz q[0],q[37];
t q[10];
t q[28];
t q[27];
t q[3];
x q[7];
cx q[24],q[25];
y q[36];
sdg q[24];
y q[15];
t q[3];
t q[39];
sdg q[6];
swap q[33],q[7];
t q[21];
cx q[37],q[24];
t q[17];
cz q[18],q[13];
t q[10];
x q[37];
sdg q[1];
swap q[16],q[2];
cz q[1],q[33];
cz q[11],q[2];
cx q[21],q[10];
x q[31];
x q[24];
swap q[4],q[7];
y q[38];
cx q[4],q[23];
cx q[10],q[9];
cz q[22],q[19];
t q[3];
cz q[36],q[20];
sdg q[11];
cz q[8],q[13];
cx q[33],q[38];
cz q[39],q[5];
z q[27];
h q[36];
y q[38];
t q[25];
z q[39];
t q[38];
cz q[25],q[19];
t q[0];
x q[8];
t q[10];
cx q[8],q[14];
s q[34];
swap q[13],q[19];
s q[39];
z q[12];
x q[0];
t q[19];
t q[21];
t q[10];
t q[8];
cx q[32],q[24];
z q[13];
t q[16];
cx q[7],q[13];
swap q[22],q[24];
y q[8];
x q[28];
cz q[8],q[2];
x q[10];
swap q[37],q[23];
swap q[29],q[16];
t q[4];
cz q[23],q[10];
t q[0];
cx q[28],q[3];
swap q[38],q[9];
t q[35];
cz q[5],q[22];
cx q[30],q[35];
t q[35];
h q[3];
t q[7];
x q[35];
cz q[12],q[7];
t q[27];
t q[33];
z q[38];
t q[17];
t q[2];
swap q[13],q[3];
cz q[35],q[5];
cx q[18],q[37];
cz q[27],q[11];
t q[16];
t q[23];
swap q[7],q[12];
t q[17];
swap q[21],q[19];
t q[22];
t q[37];
t q[12];
swap q[9],q[23];
y q[26];
t q[8];
t q[17];
t q[33];
cx q[21],q[32];
y q[0];
cx q[30],q[20];
cz q[27],q[0];
z q[9];
sdg q[24];
t q[27];
t q[13];
t q[30];
h q[31];
t q[35];
t q[30];
sdg q[8];
cx q[9],q[3];
cz q[6],q[23];
cz q[0],q[17];
t q[38];
t q[36];
t q[19];
z q[9];
cx q[27],q[2];
cx q[28],q[4];
cx q[17],q[6];
z q[19];
cx q[4],q[6];
t q[16];
cx q[35],q[17];
cz q[29],q[27];
t q[21];
t q[23];
swap q[38],q[35];
cx q[8],q[3];
t q[2];
s q[21];
t q[29];
t q[7];
x q[30];
cx q[20],q[12];
t q[10];
swap q[29],q[0];
t q[9];
swap q[3],q[13];
cz q[38],q[3];
cz q[22],q[17];
s q[17];
s q[12];
cz q[8],q[15];
sdg q[32];
t q[24];
t q[28];
t q[7];
cz q[23],q[31];
t q[9];
y q[19];
cz q[13],q[8];
t q[30];
t q[21];
t q[19];
t q[8];
t q[37];
sdg q[14];
cz q[15],q[39];
y q[29];
sdg q[34];
sdg q[23];
sdg q[12];
cz q[0],q[14];
t q[11];
t q[14];
h q[3];
t q[18];
t q[3];
y q[6];
t q[26];
s q[4];
cx q[1],q[15];
t q[18];
cx q[32],q[7];
sdg q[25];
swap q[35],q[34];
cx q[4],q[36];
y q[11];
sdg q[19];
t q[33];
y q[18];
h q[16];
t q[0];
t q[22];
swap q[16],q[20];
sdg q[10];
t q[20];
swap q[7],q[17];
t q[5];
sdg q[35];
t q[22];
z q[11];
s q[6];
x q[20];
cx q[5],q[27];
t q[4];
z q[9];
swap q[29],q[7];
y q[3];
t q[25];
swap q[26],q[14];
t q[28];
cx q[25],q[24];
t q[29];
cx q[36],q[37];
cx q[27],q[12];
x q[11];
cz q[28],q[13];
swap q[15],q[37];
s q[31];
cx q[5],q[29];
h q[24];
t q[1];
cz q[28],q[27];
t q[13];
t q[15];
cx q[13],q[33];
t q[18];
cx q[9],q[25];
t q[3];
cz q[8],q[29];
swap q[25],q[17];
x q[7];
t q[21];
z q[32];
cx q[11],q[15];
cx q[12],q[15];
sdg q[32];
cx q[17],q[6];
t q[24];
swap q[34],q[0];
swap q[29],q[31];
h q[2];
sdg q[38];
t q[5];
cz q[19],q[29];
cx q[22],q[31];
cx q[16],q[31];
cx q[29],q[12];
cz q[28],q[10];
x q[12];
h q[18];
t q[29];
s q[28];
t q[33];
z q[9];
cz q[11],q[21];
x q[3];
cx q[9],q[24];
sdg q[28];
swap q[37],q[38];
s q[14];
t q[5];
x q[24];
t q[9];
t q[36];
t q[1];
t q[33];
t q[31];
cz q[29],q[36];
cz q[31],q[33];
swap q[28],q[34];
swap q[6],q[29];
t q[31];
t q[8];
cx q[35],q[23];
t q[35];
swap q[35],q[20];
t q[4];
t q[18];
cx q[38],q[16];
t q[28];
x q[29];
t q[21];
t q[36];
t q[11];
t q[13];
t q[24];
t q[37];
swap q[11],q[0];
cz q[21],q[30];
y q[26];
t q[23];
t q[38];
x q[9];
s q[9];
s q[27];
t q[26];
cx q[32],q[1];
x q[17];
cz q[36],q[21];
t q[11];
cx q[35],q[15];
z q[26];
cz q[17],q[10];
sdg q[17];
t q[30];
t q[24];
cz q[9],q[25];
s q[11];
h q[35];
z q[26];
t q[33];
t q[10];
swap q[16],q[10];
z q[22];
t q[7];
t q[8];
cx q[27],q[19];
sdg q[12];
swap q[12],q[17];
t q[11];
cz q[28],q[27];
cx q[18],q[35];
t q[21];
sdg q[25];
cz q[27],q[12];
h q[26];
t q[19];
y q[33];
t q[31];
x q[15];
z q[9];
t q[32];
cx q[24],q[33];
cx q[7],q[25];
s q[34];
x q[18];
cx q[6],q[19];
sdg q[6];
x q[23];
cx q[25],q[6];
z q[34];
t q[27];
s q[27];
x q[11];
h q[25];
t q[31];
cz q[37],q[31];
cx q[9],q[13];
z q[4];
cx q[18],q[16];
swap q[13],q[35];
cz q[2],q[13];
t q[12];
t q[3];
t q[30];
t q[7];
sdg q[2];
swap q[32],q[11];
cx q[26],q[1];
s q[15];
sdg q[29];
t q[5];
z q[18];
cx q[33],q[7];
cz q[20],q[16];
t q[12];
cz q[9],q[38];
t q[11];
t q[38];
t q[28];
x q[19];
x q[1];
cz q[12],q[18];
z q[3];
t q[2];
t q[2];
cz q[38],q[25];
y q[24];
t q[1];
t q[19];
sdg q[4];
swap q[19],q[18];
s q[33];
t q[26];
cz q[3],q[37];
t q[20];
sdg q[33];
t q[38];
t q[22];
cz q[9],q[27];
x q[25];
t q[32];
cx q[2],q[25];
t q[1];
cx q[5],q[32];
cx q[26],q[36];
z q[20];
x q[28];
y q[16];
x q[15];
cz q[32],q[20];
t q[23];
t q[19];
swap q[12],q[4];
cz q[13],q[39];
t q[39];
cx q[25],q[15];
s q[5];
y q[18];
sdg q[23];
t q[6];
cz q[4],q[22];
sdg q[11];
t q[5];
t q[16];
s q[23];
s q[18];
t q[26];
h q[0];
t q[20];
s q[25];
cz q[22],q[12];
x q[25];
t q[17];
t q[36];
cz q[18],q[11];
cx q[37],q[25];
cx q[16],q[19];
t q[39];
cx q[8],q[7];
s q[23];
z q[18];
cz q[34],q[26];
t q[25];
cz q[26],q[28];
cz q[36],q[9];
y q[2];
y q[38];
t q[25];
sdg q[27];
h q[26];
x q[2];
t q[13];
cx q[8],q[31];
cz q[25],q[7];
x q[7];
z q[16];
sdg q[3];
z q[29];
cx q[21],q[8];
t q[2];
swap q[17],q[23];
h q[17];
t q[13];
cx q[18],q[28];
t q[4];
s q[13];
t q[35];
t q[27];
t q[35];
cz q[2],q[16];
t q[31];
s q[15];
t q[21];
t q[37];
z q[3];
sdg q[13];
t q[39];
t q[6];
x q[28];
cz q[34],q[9];
y q[16];
z q[9];
swap q[2],q[24];
x q[26];
cz q[8],q[39];
y q[33];
cz q[23],q[34];
sdg q[1];
cz q[20],q[31];
cx q[38],q[36];
x q[4];
y q[24];
t q[7];
swap q[35],q[34];
t q[39];
t q[9];
cx q[18],q[31];
t q[8];
t q[39];
cz q[34],q[39];
cx q[13],q[36];
z q[32];
s q[26];
z q[22];
cz q[22],q[17];
x q[12];
t q[33];
t q[30];
cx q[36],q[39];
t q[31];
cz q[21],q[3];
t q[1];
z q[34];
t q[33];
t q[1];
sdg q[7];
t q[15];
t q[11];
z q[1];
swap q[31],q[7];
t q[15];
z q[35];
cz q[37],q[31];
y q[19];
cz q[27],q[18];
x q[14];
sdg q[36];
x q[37];
y q[38];
t q[39];
cz q[0],q[29];
t q[7];
t q[20];
t q[11];
t q[28];
t q[6];
cx q[39],q[10];
cz q[24],q[39];
t q[2];
t q[34];
s q[27];
cx q[19],q[22];
t q[4];
t q[29];
t q[8];
t q[12];
t q[26];
t q[11];
h q[6];
t q[4];
t q[11];
cx q[34],q[30];
t q[36];
cz q[12],q[21];
y q[3];
sdg q[0];
t q[26];
t q[6];
t q[35];
t q[12];
cx q[12],q[31];
swap q[12],q[7];
y q[24];
cz q[28],q[36];
cx q[2],q[37];
swap q[13],q[34];
t q[39];
x q[7];
s q[39];
h q[35];
cx q[9],q[26];
t q[22];
x q[22];
cz q[15],q[35];
cz q[30],q[27];
x q[5];
cx q[5],q[29];
s q[36];
z q[7];
y q[19];
cz q[0],q[19];
swap q[14],q[13];
cz q[9],q[39];
cx q[16],q[1];
t q[2];
t q[18];
t q[16];
cz q[18],q[31];
t q[16];
h q[33];
t q[17];
t q[4];
swap q[7],q[25];
swap q[2],q[22];
t q[35];
cz q[10],q[19];
cx q[3],q[11];
t q[9];
cz q[34],q[21];
x q[3];
cz q[37],q[18];
y q[21];
swap q[11],q[17];
s q[15];
t q[22];
t q[12];
s q[28];
swap q[9],q[22];
t q[27];
t q[30];
z q[35];
cx q[39],q[17];
x q[19];
x q[31];
y q[8];
s q[26];
cx q[21],q[5];
swap q[28],q[9];
cz q[28],q[24];
cz q[16],q[31];
t q[14];
swap q[24],q[21];
t q[6];
t q[22];
swap q[32],q[15];
t q[1];
t q[38];
cz q[24],q[22];
swap q[32],q[8];
t q[16];
t q[12];
swap q[35],q[8];
t q[38];
x q[39];
cz q[8],q[36];
t q[20];
h q[30];
t q[19];
t q[0];
s q[36];
swap q[15],q[23];
y q[32];
t q[30];
cz q[22],q[7];
t q[30];
z q[1];
s q[30];
sdg q[31];
t q[1];
t q[8];
y q[22];
t q[22];
cz q[12],q[26];
y q[17];
t q[0];
cz q[28],q[16];
x q[23];
swap q[13],q[36];
x q[11];
t q[4];
t q[17];
t q[33];
x q[8];
cz q[19],q[26];
z q[26];
swap q[9],q[37];
cz q[27],q[29];
y q[6];
t q[12];
t q[1];
swap q[9],q[38];
s q[25];
s q[6];
y q[10];
s q[11];
x q[34];
t q[8];
swap q[3],q[8];
swap q[9],q[0];
t q[29];
h q[3];
h q[16];
x q[13];
t q[32];
t q[1];
cz q[38],q[34];
cz q[15],q[22];
t q[4];
swap q[8],q[30];
swap q[8],q[18];
z q[27];
t q[1];
t q[39];
s q[36];
t q[4];
x q[7];
t q[29];
t q[25];
cz q[22],q[7];
x q[4];
y q[2];
x q[19];
t q[38];
y q[37];
s q[13];
h q[3];
t q[10];
cx q[37],q[24];
swap q[13],q[12];
y q[11];
t q[29];
t q[37];
swap q[26],q[7];
swap q[23],q[16];
t q[13];
swap q[13],q[29];
y q[32];
t q[30];
swap q[5],q[20];
cx q[37],q[16];
t q[38];
h q[26];
t q[7];
z q[31];
t q[29];
x q[10];
t q[27];
cz q[17],q[35];
cz q[2],q[10];
t q[4];
t q[16];
t q[6];
swap q[18],q[31];
x q[32];
sdg q[21];
t q[35];
h q[34];
swap q[39],q[32];
cz q[39],q[4];
y q[10];
t q[9];
t q[16];
t q[11];
cx q[33],q[7];
sdg q[27];
cx q[3],q[5];
t q[1];
y q[7];
z q[16];
x q[1];
swap q[12],q[4];
z q[5];
cz q[31],q[8];
cx q[9],q[28];
h q[3];
cx q[33],q[38];
t q[5];
s q[24];
cz q[11],q[18];
cz q[21],q[38];
t q[18];
h q[28];
cx q[14],q[39];
t q[39];
t q[39];
t q[19];
t q[29];
t q[24];
t q[1];
t q[36];
y q[2];
t q[12];
t q[28];
t q[2];
t q[4];
t q[19];
t q[16];
swap q[36],q[9];
cz q[5],q[0];
y q[28];
cz q[29],q[15];
t q[16];
t q[25];
y q[22];
t q[22];
t q[33];
cx q[34],q[31];
cx q[0],q[31];
y q[10];
s q[26];
t q[22];
t q[15];
y q[33];
x q[24];
t q[16];
t q[18];
swap q[37],q[36];
h q[21];
sdg q[16];
sdg q[37];
cx q[14],q[7];
t q[20];
t q[27];
swap q[35],q[0];
cz q[13],q[7];
cx q[30],q[14];
x q[32];
swap q[4],q[34];
h q[9];
y q[34];
sdg q[6];
sdg q[6];
cx q[36],q[15];
t q[34];
y q[22];
t q[33];
t q[38];
t q[22];
t q[22];
swap q[18],q[36];
swap q[8],q[12];
sdg q[25];
t q[33];
t q[31];
s q[32];
t q[22];
z q[37];
t q[37];
t q[27];
cz q[26],q[25];
h q[20];
t q[24];
sdg q[32];
t q[25];
cx q[26],q[31];
x q[12];
z q[30];
t q[37];
t q[7];
x q[6];
sdg q[28];
x q[34];
t q[39];
swap q[4],q[5];
cx q[14],q[1];
cz q[2],q[33];
cx q[34],q[3];
t q[35];
t q[16];
swap q[17],q[14];
cx q[1],q[3];
t q[23];
t q[38];
sdg q[20];
h q[23];
cz q[1],q[9];
cx q[34],q[3];
t q[24];
z q[24];
swap q[0],q[9];
t q[38];
sdg q[32];
swap q[9],q[18];
cx q[25],q[38];
cz q[36],q[15];
t q[34];
h q[0];
x q[7];
t q[4];
t q[23];
swap q[1],q[7];
t q[10];
s q[34];
sdg q[2];
t q[38];
t q[25];
h q[20];
cx q[14],q[25];
h q[27];
cx q[13],q[5];
swap q[37],q[35];
cx q[16],q[28];
sdg q[7];
cx q[21],q[39];
t q[0];
cz q[36],q[7];
cx q[22],q[24];
cz q[2],q[7];
swap q[29],q[32];
cx q[28],q[31];
y q[13];
z q[30];
t q[28];
t q[16];
t q[10];
t q[36];
t q[12];
y q[36];
t q[12];
t q[6];
t q[26];
t q[32];
t q[35];
x q[35];
cz q[11],q[21];
cz q[28],q[6];
cx q[16],q[34];
t q[22];
swap q[38],q[20];
t q[31];
cx q[24],q[5];
cz q[9],q[36];
sdg q[33];
s q[25];
t q[18];
y q[2];
cx q[34],q[12];