%CNV detection from exon sequencing data
%model min_(v,d,C) 1/2*||Y/v-dC||^2+la*||C||_0
%junbo.duan@mail.xjtu.edu.cn
%2014/11/15

close all;clear all;clc;

%%%%%%%%%%%%%%%%%%parameter configuration%%%%%%
p=0.92;
len_min=3;
gap_min=4;
num_it=100;

%%%%%%%%%%%%%%%%%data loading%%%%%%%%%%%%%%%%%
fid=fopen('G:\linux_suse\share_linux_suse\1000exon\sampledata\probes.txt');
probes=textscan(fid,'%n %n %n %s');%loading probes
fclose(fid);

tab_files=getAllFiles('G:\linux_suse\share_linux_suse\1000exon\sampledata\RPKM_data\');

for i=1:length(tab_files)
    fid=fopen(tab_files{i});
    datatab=textscan(fid,'%f %f %f');%loading data from files
    fclose(fid);
    rd(:,i)=datatab{2};
    rpkm(:,i)=datatab{3};
end

idx_sel=(median(rpkm,2)>1)&(probes{1}<=22);%quality control and autosomes

Y=rpkm(idx_sel,:);
%Y=rd(idx_sel,:);
[M,N]=size(Y);

%%%%%%%%%%%%%%%%%%%%%%optimization%%%%%%%%%%%%%%%%%%%%%%%
for n=1:N
    v(n)=norm(Y(:,n));%initializing v
end
v=v(:);

for k=1:num_it%starting optimization
    
    k
    tic;
    
    tab_d=zeros(1,N);%exhaustive search to solve 1D Lp optimization
    for ii=1:M
        tab_d=sum((abs(ones(N,1)*(Y(ii,:)./(v'))-Y(ii,:)'*ones(1,N))).^p,2);
        [NA,idx]=min(tab_d);
        d(ii)=Y(ii,idx)/v(idx);
    end
    
    d=d(:)/norm(d)*(mean(std(Y')))*M;
    
    for n=1:N
        w=(abs(d-Y(:,n)*v(n)^(-1))).^(p-2);
        w(w==inf)=0;%set the .^(p-2) of zero elements to zero
        v(n)=sum(w.*Y(:,n).*Y(:,n))*(sum(w.*Y(:,n).*d))^(-1);
    end
    
    toc;
end


%%%%%%%%%%%%%%%%%%%%%C estimation%%%%%%%%%%%%%%%%
rsd=zeros(size(Y));
C=zeros(size(Y));

for i=1:N
    rsd(:,i)=(abs(Y(:,i)/v(i)-d)).^p;
    C(:,i)=(Y(:,i)/v(i))./d;
end

mu=(mean(std((Y*inv(diag(v)))')))^p*log(M*N)*1;%penalty value estimation

idx_C=find(rsd<mu);
C(idx_C)=1;

CO=NaN*ones(size(rpkm));%C is the masked verson of CO, mask is idx_sel.
CO(idx_sel,:)=C;

%%%%%%%%%%%%%%%%%%%CNV calling%%%%%%%%%%%%%%%%%%%
for i=1:N
    CNV{i}=[];
    idx=find(CO(:,i)~=1);
    
    j=1;
    while j<length(idx)
        len=1;
        while (sign(CO(idx(j+len),i))==sign(CO(idx(j),i)))...
                &&((idx(j+len)-idx(j))==len)...
                &&((j+len)<length(idx))...
                &&(probes{1}(idx(j))==probes{1}(idx(j+len)))
            len=len+1;
        end
        if len>=len_min
            CNV{i}=[CNV{i};[idx(j),idx(j+len-1)]];
        end
        j=j+len;
    end
    
    for j=size(CNV{i},1)-1:-1:1
        if (CNV{i}(j+1,1)-CNV{i}(j,2))<=gap_min
            CNV{i}(j,2)=CNV{i}(j+1,2);
            CNV{i}(j+1,:)=[];
        end
    end
    
    for j=1:size(CNV{i})
        CNV{i}(j,3)=mean(CO(CNV{i}(j,1):CNV{i}(j,2),i))*2;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%output saving%%%%%%%%%%%%%%%%%
fid=fopen('G:\linux_suse\share_linux_suse\1000exon\outputs.txt','w');
fprintf(fid,'CNV\tsample\tchromosome\tstart\tend\tcopy estimate\n');

k=1;
for i=1:N
    str_file=tab_files{i};
    loc=strfind(str_file,'\NA1');
    id_sample=str2num(str_file(loc+3:loc+7));
    for j=1:size(CNV{i},1)
        chr=probes{1}(CNV{i}(j,1));
        bpstart=probes{2}(CNV{i}(j,1));
        bpend=probes{3}(CNV{i}(j,2));
        cp=CNV{i}(j,3);
        
        fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%f\n',k,id_sample,chr,bpstart,bpend,cp);
        k=k+1;
    end
end

fclose(fid);