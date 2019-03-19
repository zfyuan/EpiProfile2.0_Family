function mergeC12C13
%%

clc;

file1 = 'histone_ratios.xls';
file2 = 'histone_ratios_C13.xls';
file3 = 'histone_ratios_C12+C13.xls';

strs1 = read_file(file1);
strs2 = read_file(file2);

[lineno1,area1,rt1] = get_area_rt(strs1);
[lineno2,area2] = get_area(strs2);

area3 = merge_area_C12C13(lineno1,area1,lineno2,area2);

ratio3 = get_ratio(area3,lineno1);

output(file3,strs1,lineno1,ratio3,area3,rt1);

function strs1 = read_file(file1)
%%

strs1 = {};
no = 0;
fp = fopen(file1,'r');
if -1==fp
    fprintf(1,'can not open: %s\n',file1);
    return;
end

while 1
    str = fgetl(fp);
    if 1==isempty(str)
        continue;
    end
    
    no = no + 1;
    strs1{no,1} = str;%#ok
    
    if 1==feof(fp)
        break;
    end
end

fclose(fp);

function [lineno1,area1,rt1] = get_area_rt(strs1)
%%

p = strfind(strs1{1,1},',');
nfiles = length(p)/3;

lineno1 = [];
pnum = 0;
for ino=4:size(strs1,1)
    c_str = strs1{ino,1};
    if 1==strcmp(c_str(1:2),'H2') || 1==strcmp(c_str(1:2),'H3') || 1==strcmp(c_str(1:2),'H4')
        pnum = pnum + 1;
        lineno1(pnum,1) = ino;%#ok
    end
end

area1 = zeros([pnum,nfiles]);
rt1 = zeros([pnum,nfiles]);
key = '	';% \t

for pno=1:pnum
    ino = lineno1(pno,1);
    c_str = strs1{ino,1};
    pp = [strfind(c_str,key),length(c_str)+1];
    for jno=nfiles+2:nfiles*2+1
        qno = jno-(nfiles+1);
        area1(pno,qno) = str2num(c_str(pp(jno)+1:pp(jno+1)-1));%#ok
    end
    for jno=nfiles*2+3:nfiles*3+2
        qno = jno-(nfiles*2+2);
        rt1(pno,qno) = str2num(c_str(pp(jno)+1:pp(jno+1)-1));%#ok
    end
end

function [lineno2,area2] = get_area(strs2)
%%

p = strfind(strs2{1,1},',');
nfiles = length(p)/3;

lineno2 = [];
pnum = 0;
for ino=4:size(strs2,1)
    c_str = strs2{ino,1};
    if 0==isempty(strfind(c_str,'C2'))%#ok
        pnum = pnum + 1;
        lineno2(pnum,1) = ino;%#ok
    end
end

area2 = zeros([pnum,nfiles]);
key = '	';% \t

for pno=1:pnum
    ino = lineno2(pno,1);
    c_str = strs2{ino,1};
    pp = [strfind(c_str,key),length(c_str)+1];
    for jno=nfiles+2:nfiles*2+1
        qno = jno-(nfiles+1);
        area2(pno,qno) = str2num(c_str(pp(jno)+1:pp(jno+1)-1));%#ok
    end
end

function area3 = merge_area_C12C13(lineno1,area1,lineno2,area2)
%%

area3 = area1;
nfiles = size(area1,2);

% H3_9_17
x1 = find(lineno1==8);
x2 = find(lineno1==9);
xx1 = find(lineno2==6);
xx2 = find(lineno2==7);
a1 = area1(x1,:);
a2 = area1(x2,:);
sum1 = a1+a2+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==10);
xx1 = find(lineno2==10);
area3(x1,:) = area3(x1,:) + area2(xx1,:);%#ok

x1 = find(lineno1==11);
xx1 = find(lineno2==13);
area3(x1,:) = area3(x1,:) + area2(xx1,:);%#ok

x1 = find(lineno1==12);
xx1 = find(lineno2==16);
area3(x1,:) = area3(x1,:) + area2(xx1,:);%#ok

x1 = find(lineno1==13);
xx1 = find(lineno2==19);
xx2 = find(lineno2==20);
xx3 = find(lineno2==21);
area3(x1,:) = area3(x1,:) + area2(xx1,:) + area2(xx2,:) + area2(xx3,:);%#ok

% H3_18_26
x1 = find(lineno1==19);
x2 = find(lineno1==20);
xx1 = find(lineno2==25);
xx2 = find(lineno2==26);
a1 = area1(x1,:);
a2 = area1(x2,:);
sum1 = a1+a2+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==21);
xx1 = find(lineno2==29);
xx2 = find(lineno2==30);
xx3 = find(lineno2==31);
area3(x1,:) = area3(x1,:) + area2(xx1,:) + area2(xx2,:) + area2(xx3,:);%#ok

% H3_27_40
x1 = find(lineno1==37);
xx1 = find(lineno2==34);
area3(x1,:) = area3(x1,:) + area2(xx1,:);%#ok

% H33_27_40
x1 = find(lineno1==53);
xx1 = find(lineno2==37);
area3(x1,:) = area3(x1,:) + area2(xx1,:);%#ok

% H4_4_17
x1 = find(lineno1==56);
x2 = find(lineno1==57);
x3 = find(lineno1==58);
x4 = find(lineno1==59);
xx1 = find(lineno2==43);
xx2 = find(lineno2==44);
xx3 = find(lineno2==45);
xx4 = find(lineno2==46);
a1 = area1(x1,:);
a2 = area1(x2,:);
a3 = area1(x3,:);
a4 = area1(x4,:);
sum1 = a1+a2+a3+a4+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:)+area2(xx3,:)+area2(xx4,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
    area3(x3,jno) = area3(x3,jno) + a3(jno)/sum1(jno)*sum2(jno);
    area3(x4,jno) = area3(x4,jno) + a4(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==60);
x2 = find(lineno1==61);
x3 = find(lineno1==62);
x4 = find(lineno1==63);
x5 = find(lineno1==64);
x6 = find(lineno1==65);
xx1 = find(lineno2==49);
xx2 = find(lineno2==50);
a1 = area1(x1,:);
a2 = area1(x2,:);
a3 = area1(x3,:);
a4 = area1(x4,:);
a5 = area1(x5,:);
a6 = area1(x6,:);
sum1 = a1+a2+a3+a4+a5+a6+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
    area3(x3,jno) = area3(x3,jno) + a3(jno)/sum1(jno)*sum2(jno);
    area3(x4,jno) = area3(x4,jno) + a4(jno)/sum1(jno)*sum2(jno);
    area3(x5,jno) = area3(x5,jno) + a5(jno)/sum1(jno)*sum2(jno);
    area3(x6,jno) = area3(x6,jno) + a6(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==66);
x2 = find(lineno1==67);
x3 = find(lineno1==68);
x4 = find(lineno1==69);
xx1 = find(lineno2==53);
xx2 = find(lineno2==54);
xx3 = find(lineno2==55);
a1 = area1(x1,:);
a2 = area1(x2,:);
a3 = area1(x3,:);
a4 = area1(x4,:);
sum1 = a1+a2+a3+a4+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:)+area2(xx3,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
    area3(x3,jno) = area3(x3,jno) + a3(jno)/sum1(jno)*sum2(jno);
    area3(x4,jno) = area3(x4,jno) + a4(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==70);
xx1 = find(lineno2==58);
xx2 = find(lineno2==59);
xx3 = find(lineno2==60);
xx4 = find(lineno2==61);
area3(x1,:) = area3(x1,:) + area2(xx1,:) + area2(xx2,:) + area2(xx3,:) + area2(xx4,:);%#ok

% H2A1_4_11
x1 = find(lineno1==73);
x2 = find(lineno1==74);
xx1 = find(lineno2==65);
xx2 = find(lineno2==66);
a1 = area1(x1,:);
a2 = area1(x2,:);
sum1 = a1+a2+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==75);
xx1 = find(lineno2==69);
xx2 = find(lineno2==70);
xx3 = find(lineno2==71);
area3(x1,:) = area3(x1,:) + area2(xx1,:) + area2(xx2,:) + area2(xx3,:);%#ok

% H2AV_1_19
x1 = find(lineno1==80);
x2 = find(lineno1==81);
x3 = find(lineno1==82);
x4 = find(lineno1==83);
xx1 = find(lineno2==77);
xx2 = find(lineno2==78);
xx3 = find(lineno2==79);
xx4 = find(lineno2==80);
a1 = area1(x1,:);
a2 = area1(x2,:);
a3 = area1(x3,:);
a4 = area1(x4,:);
sum1 = a1+a2+a3+a4+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:)+area2(xx3,:)+area2(xx4,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
    area3(x3,jno) = area3(x3,jno) + a3(jno)/sum1(jno)*sum2(jno);
    area3(x4,jno) = area3(x4,jno) + a4(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==84);
x2 = find(lineno1==85);
x3 = find(lineno1==86);
x4 = find(lineno1==87);
x5 = find(lineno1==88);
x6 = find(lineno1==89);
xx1 = find(lineno2==83);
xx2 = find(lineno2==84);
a1 = area1(x1,:);
a2 = area1(x2,:);
a3 = area1(x3,:);
a4 = area1(x4,:);
a5 = area1(x5,:);
a6 = area1(x6,:);
sum1 = a1+a2+a3+a4+a5+a6+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
    area3(x3,jno) = area3(x3,jno) + a3(jno)/sum1(jno)*sum2(jno);
    area3(x4,jno) = area3(x4,jno) + a4(jno)/sum1(jno)*sum2(jno);
    area3(x5,jno) = area3(x5,jno) + a5(jno)/sum1(jno)*sum2(jno);
    area3(x6,jno) = area3(x6,jno) + a6(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==90);
x2 = find(lineno1==91);
x3 = find(lineno1==92);
x4 = find(lineno1==93);
xx1 = find(lineno2==87);
xx2 = find(lineno2==88);
xx3 = find(lineno2==89);
a1 = area1(x1,:);
a2 = area1(x2,:);
a3 = area1(x3,:);
a4 = area1(x4,:);
sum1 = a1+a2+a3+a4+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:)+area2(xx3,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
    area3(x3,jno) = area3(x3,jno) + a3(jno)/sum1(jno)*sum2(jno);
    area3(x4,jno) = area3(x4,jno) + a4(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==94);
xx1 = find(lineno2==92);
xx2 = find(lineno2==93);
xx3 = find(lineno2==94);
xx4 = find(lineno2==95);
area3(x1,:) = area3(x1,:) + area2(xx1,:) + area2(xx2,:) + area2(xx3,:) + area2(xx4,:);%#ok

% H2AZ_1_19
x1 = find(lineno1==97);
x2 = find(lineno1==98);
x3 = find(lineno1==99);
x4 = find(lineno1==100);
xx1 = find(lineno2==101);
xx2 = find(lineno2==102);
xx3 = find(lineno2==103);
xx4 = find(lineno2==104);
a1 = area1(x1,:);
a2 = area1(x2,:);
a3 = area1(x3,:);
a4 = area1(x4,:);
sum1 = a1+a2+a3+a4+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:)+area2(xx3,:)+area2(xx4,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
    area3(x3,jno) = area3(x3,jno) + a3(jno)/sum1(jno)*sum2(jno);
    area3(x4,jno) = area3(x4,jno) + a4(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==101);
x2 = find(lineno1==102);
x3 = find(lineno1==103);
x4 = find(lineno1==104);
x5 = find(lineno1==105);
x6 = find(lineno1==106);
xx1 = find(lineno2==107);
xx2 = find(lineno2==108);
a1 = area1(x1,:);
a2 = area1(x2,:);
a3 = area1(x3,:);
a4 = area1(x4,:);
a5 = area1(x5,:);
a6 = area1(x6,:);
sum1 = a1+a2+a3+a4+a5+a6+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
    area3(x3,jno) = area3(x3,jno) + a3(jno)/sum1(jno)*sum2(jno);
    area3(x4,jno) = area3(x4,jno) + a4(jno)/sum1(jno)*sum2(jno);
    area3(x5,jno) = area3(x5,jno) + a5(jno)/sum1(jno)*sum2(jno);
    area3(x6,jno) = area3(x6,jno) + a6(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==107);
x2 = find(lineno1==108);
x3 = find(lineno1==109);
x4 = find(lineno1==110);
xx1 = find(lineno2==111);
xx2 = find(lineno2==112);
xx3 = find(lineno2==113);
a1 = area1(x1,:);
a2 = area1(x2,:);
a3 = area1(x3,:);
a4 = area1(x4,:);
sum1 = a1+a2+a3+a4+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:)+area2(xx3,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
    area3(x3,jno) = area3(x3,jno) + a3(jno)/sum1(jno)*sum2(jno);
    area3(x4,jno) = area3(x4,jno) + a4(jno)/sum1(jno)*sum2(jno);
end

x1 = find(lineno1==111);
xx1 = find(lineno2==116);
xx2 = find(lineno2==117);
xx3 = find(lineno2==118);
xx4 = find(lineno2==119);
area3(x1,:) = area3(x1,:) + area2(xx1,:) + area2(xx2,:) + area2(xx3,:) + area2(xx4,:);%#ok

% H2A1_12_17
x1 = find(lineno1==114);
x2 = find(lineno1==115);
xx1 = find(lineno2==123);
xx2 = find(lineno2==124);
a1 = area1(x1,:);
a2 = area1(x2,:);
sum1 = a1+a2+repmat(eps,[1,nfiles]);
sum2 = area2(xx1,:)+area2(xx2,:);%#ok
for jno=1:nfiles
    area3(x1,jno) = area3(x1,jno) + a1(jno)/sum1(jno)*sum2(jno);
    area3(x2,jno) = area3(x2,jno) + a2(jno)/sum1(jno)*sum2(jno);
end

function ratio3 = get_ratio(area3,lineno1)
%%

a = lineno1(2:end);
b = lineno1(1:end-1);

c = a-b;
c(end+1) = 2;

idx = [0;find(c==2)];

[pnum,nfiles] = size(area3);
ratio3 = zeros([pnum,nfiles]);

for ino=1:length(idx)-1
    IX = idx(ino)+1:idx(ino+1);
    for jno=1:nfiles
        sum1 = sum(area3(IX,jno))+eps;
        ratio3(IX,jno) = area3(IX,jno)/sum1;
    end
end

function output(file3,strs1,lineno1,ratio3,area3,rt1)
%%

fp = fopen(file3,'w');
if -1==fp
    fprintf(1,'can not open: %s\n',file3);
    return;
end

nfiles = size(ratio3,2);
key = '	';% \t
for ino=1:size(strs1,1)
    if 0==ismember(ino,lineno1)
        fprintf(fp,'%s\n',strs1{ino,1});
    else
        c_str = strs1{ino,1};
        p = strfind(c_str,key);
        prefix = c_str(1:p(1)-1);
        fprintf(fp,'%s',prefix);
        x1 = find(lineno1==ino);
        for jno=1:nfiles
            fprintf(fp,'\t%f',ratio3(x1,jno));
        end
        fprintf(fp,'\t');
        for jno=1:nfiles
            fprintf(fp,'\t%e',area3(x1,jno));
        end
        fprintf(fp,'\t');
        for jno=1:nfiles
            fprintf(fp,'\t%.2f',rt1(x1,jno));
        end
        fprintf(fp,'\n');
    end
end

fclose(fp);