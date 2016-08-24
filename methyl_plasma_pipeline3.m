function methyl_plasma_pipeline3()
%% Pipeline for plasma methylation analysis
%% 
%% The user will be prompted for three files: a sequence fastq file, a barcode fastq file, and an input tab-delimited text file.
%% The input text file should have the following format: target_name sample_name barcode pre-bis_amplicon
%% Output: fastq files per barcode, raw,histogram,summary file for each target-sample, smmary of all results to Summary_results.txt 
%% For questions please contact joshua.moss@mail.huji.ac.il

%matlabpool open local
num_cores = feature('numCores');
%num_cores = 32;
parpool('local',num_cores);

%% Files enter
%%start_folder = '/mnt/lustre/hms-01/fs01'
%[indexFile,folder,~] = uigetfile('./*.fastq','Choose Index file');
%[readsFile,~,~] = uigetfile([folder '/*.fastq'],'Choose Reads file');
%[seqFile,~,~] = uigetfile([folder '/*.txt'],'Choose Sequence file');

folder = './';
seqFile = 'sample_sheet.txt';
readsFile = 'Undetermined_S0_R1_001.fastq';
indexFile = 'Undetermined_S0_I1_001.fastq';

cd(folder);

%% import info from sequence file

fid = fopen(seqFile);
info = textscan(fid,'%s %s %s %s','Delimiter',{'\t'},'CollectOutput',1,'Whitespace',' \b"');
fclose(fid);
info = info{1};
sample_name1 = info(:,1);
sample_name2 = info(:,2);
sample_barcodes = info(:,3);
sample_seq = info(:,4);
[barcodes,~,sample_barcodes_num] = unique(sample_barcodes);
sampleNum = length(sample_barcodes_num);

%% write input summary to file
fid = fopen('input_summary.txt','w');
fprintf(fid, '%s\t', 'Sample #');
fprintf(fid, '%s\t', 'Name');
fprintf(fid, '%s\t', 'Info');
fprintf(fid, '%s\t', 'Barcode');
fprintf(fid, '%s\t', 'Sequence');
for i = 1:sampleNum
    fprintf(fid,'\n');
    fprintf(fid,'%f\t',i);
    for j = 1:4
        fprintf(fid,'%s\t',info{i,j});
    end
end
fclose(fid);

%% write fasta files
for i = 1:sampleNum
    fastawrite(['sample' int2str(i) '.fasta'], [sample_name1{i} ' ' sample_name2{i}], sample_seq{i});
end

%% write barcode summary
fid = fopen('barcode_order.txt','w');
fprintf(fid,'%s\t', 'Barcode #');
fprintf(fid, '%s\t', 'Barcode');
for i = 1:length(barcodes)
    fprintf(fid,'\n');
    fprintf(fid,'%f\t', i);
    fprintf(fid,'%s\t',barcodes{i});
end
fclose(fid);


%% quality filter 
%% change to allow parallel barcode division
[~, result] = system( ['wc -l ', indexFile] );
line_count = strsplit(result,' ');
numlines_I = str2num(line_count{1});

[~, result] = system( ['wc -l ', readsFile] );
line_count = strsplit(result,' ');
numlines_R = str2num(line_count{1});

if ~(numlines_I==numlines_R)
    error('Index and read fastq files are of different lengths');
end

smallfile_lines = ceil((numlines_I/num_cores)/4)*4;
status = system(['split -l ',int2str(smallfile_lines),' -d ', readsFile,' ',readsFile,'_']);
status = system(['split -l ',int2str(smallfile_lines),' -d ', indexFile,' ',indexFile,'_']);

fprintf('%s','Fasta files written. Quality filter and barcode divison will now begin. Please be patient!');

parfor i = 0:(num_cores-1)
    filter_and_seperate_par(indexFile,readsFile,barcodes,i);
end
parfor i = 1:length(barcodes)
    status = system( ['cat barcode' int2str(i) '.fastq_* > barcode' int2str(i) '.fastq'] );
end
status = system('rm -f *fastq_*');
fprintf('%s','Barcode division is done! Analysis will now be performed');
fprintf('\n');

%% Methanalysis run

% Assign Params structs for each sample
Params_array = cell(sampleNum,1);
for i = 1:sampleNum
    Params.RefFile = ['sample' int2str(i) '.fasta'];
    Params.ReadsFile = ['barcode' int2str(sample_barcodes_num(i)) '.fastq'];
    Params_array{i} = Params;
end


% RUN
parfor i = 1:sampleNum
    RunMethAnalysisFile(Params_array{i});
end

poolobj = gcp('nocreate');
delete(poolobj);

%% Summary file
num_CpGs = zeros(sampleNum,1);
for i = 1:sampleNum
    num_CpGs(i) = length(strfind(sample_seq{i},'CG'));
end
max_num_CpGs = max(num_CpGs);

numT = zeros(sampleNum,max_num_CpGs+3);
numT(:,:) = NaN;
for i = 1:sampleNum
    numT_sample = numToutput(num_CpGs(i),['sample' int2str(i)]);
    numT(i,1:length(numT_sample)) = numT_sample;
end

%print
fid = fopen('Summary_results.txt','w');
fprintf(fid,'%s\t', 'Sample #');
fprintf(fid,'%s\t', 'Gene');
fprintf(fid,'%s\t', 'Sample');
fprintf(fid,'%s\t', 'CpGs');
fprintf(fid,'%s\t', 'All reads');
fprintf(fid,'%s\t', 'All T');
for i = 1:max_num_CpGs
    fprintf(fid,'%s\t', ['All T - ' int2str(i)]);
end

for i = 1:sampleNum
    fprintf(fid,'\n');
    fprintf(fid,'%f\t', i);
    fprintf(fid,'%s\t', sample_name1{i});
    fprintf(fid,'%s\t', sample_name2{i});    
    for j = 1:length(numT(1,:))
        if isnan(numT(i,j))
            break;
        end
        fprintf(fid, '%f\t', numT(i,j));
    end
end
fclose(fid);  
%msgbox('Congratulations - All done!');
end
function filter_and_seperate(indexFile,readsFile,barcodes)
%% Filters reads by the quality of their barcode read and seperates reads by barcode

% set quality cutoff (can be adjusted)
qt = 32;

% read in relevant barcodes
for i = 1:length(barcodes)
    barcodes{i} = seqrcomplement(barcodes{i});
end

% get general sequence information
info_I = fastqinfo(indexFile);
len_I = double(info_I.NumberOfEntries);

fid_barcodes = cell(length(barcodes),1);
for i = 1:length(barcodes)
    fid_barcodes{i} = fopen(['barcode' int2str(i) '.fastq'],'a');
end

fidI = fopen(indexFile);
fidR = fopen(readsFile);

tmp = fgetl(fidI);
i=0;
tic
while tmp ~= -1
    i = i+1;
    if (mod(i,1000)==0)    
        t = toc;
        rem = t*(len_I-i)/i;    
        disp(['Seperating barcodes ' num2str(i) '/' num2str(len_I) '(' num2str(round(100*i/len_I)) '%) Remaining: ' datestr(datenum(0,0,0,0,0,rem),'HH:MM:SS')]);
    end
    
    % read barcode info
    seq_I = fgetl(fidI);
    tmp = fgetl(fidI);
    qual_I = fgetl(fidI);
    tmp = fgetl(fidI);
    bc_compare = find(strcmp(seq_I,barcodes));
    
    % read seq info
    header_R = fgetl(fidR);
    seq_R = fgetl(fidR);
    header2_R = fgetl(fidR);
    qual_R = fgetl(fidR);    
    
    if ~isempty(bc_compare)
        % filter quality barcode reads
        qual_I_num = mean(double(qual_I)-33);
        if qual_I_num > qt
            % write to barcode file
            %fidOut = ;
            fprintf(fid_barcodes{bc_compare}, '%s\n', header_R);
            fprintf(fid_barcodes{bc_compare}, '%s\n', seq_R);
            fprintf(fid_barcodes{bc_compare}, '%s\n', header2_R);
            fprintf(fid_barcodes{bc_compare}, '%s\n', qual_R);
        end        
    end    
      
end

fclose(fidI);
fclose(fidR);

for i = 1:length(barcodes)
    fclose(fid_barcodes{i});
end

fprintf('%s','Quality filter is done.');
fprintf('\n');

end
function filter_and_seperate_par(indexFile,readsFile,barcodes,cur_i)
%% Filters reads by the quality of their barcode read and seperates reads by barcode

if cur_i<10
    cur_i = ['0' int2str(cur_i)];
else
    cur_i = int2str(cur_i);
end
indexFile = [indexFile '_' cur_i];
readsFile = [readsFile '_' cur_i];
    
% set quality cutoff (can be adjusted)
qt = 32;

% read in relevant barcodes
for i = 1:length(barcodes)
    barcodes{i} = seqrcomplement(barcodes{i});
end

%% get general sequence information
%info_I = fastqinfo(indexFile);
%len_I = double(info_I.NumberOfEntries);

fid_barcodes = cell(length(barcodes),1);
for i = 1:length(barcodes)
    fid_barcodes{i} = fopen(['barcode' int2str(i) '.fastq_' cur_i],'a');
end

fidI = fopen(indexFile);
fidR = fopen(readsFile);

tmp = fgetl(fidI);
i=0;
tic
while tmp ~= -1
    i = i+1;
    %if (mod(i,1000)==0)    
    %    t = toc;
    %    rem = t*(len_I-i)/i;    
    %    disp(['Seperating barcodes ' num2str(i) '/' num2str(len_I) '(' num2str(round(100*i/len_I)) '%) Remaining: ' datestr(datenum(0,0,0,0,0,rem),'HH:MM:SS')]);
    %end
    
    % read barcode info
    seq_I = fgetl(fidI);
    tmp = fgetl(fidI);
    qual_I = fgetl(fidI);
    tmp = fgetl(fidI);
    bc_compare = find(strcmp(seq_I,barcodes));
    
    % read seq info
    header_R = fgetl(fidR);
    seq_R = fgetl(fidR);
    header2_R = fgetl(fidR);
    qual_R = fgetl(fidR);    
    
    if ~isempty(bc_compare)
        % filter quality barcode reads
        qual_I_num = mean(double(qual_I)-33);
        if qual_I_num > qt
            % write to barcode file
            %fidOut = ;
            fprintf(fid_barcodes{bc_compare}, '%s\n', header_R);
            fprintf(fid_barcodes{bc_compare}, '%s\n', seq_R);
            fprintf(fid_barcodes{bc_compare}, '%s\n', header2_R);
            fprintf(fid_barcodes{bc_compare}, '%s\n', qual_R);
        end        
    end    
      
end

fclose(fidI);
fclose(fidR);

for i = 1:length(barcodes)
    fclose(fid_barcodes{i});
end

fprintf('%s','Quality filter is done.');
fprintf('\n');

end
function RunMethAnalysisFile(Params)
%% Performs local alignment of all reads in a fastq file to target amplicon in fasta file
%% Outputs three files for each sample: 
%% .meth_raw - shows raw alignments and CpG information
%% .meth_CpG_hist - shows counts of each different combinations of CpGs in aligned reads
%% .meth_summary - shows the methylation information for each CpG from target amplicon
Params.OutputFileName = strrep(Params.RefFile,'.fasta','.meth_raw');
name_core = strrep(Params.RefFile,'.fasta','');


Params.Thresholds.PercentCorrectAlignment = 0.8; % sets minimum percent of target sequence aligned to read; can be adjusted
Params.Thresholds.MinReadLength = 50; % sets minimum read length; can be adjusted
Params.Thresholds.MinAlignmentLength = 50; % sets minimum length of aligned region in read; can be adjusted


if (exist(Params.RefFile) == 0)
    disp(['Reference file does not exist: ' Params.RefFile]);
    return;
end
if (exist(Params.ReadsFile) == 0)
    error(['Reads file does not exist: ' Params.ReadsFile]);
end

disp(['Counting Lines for file: ' Params.OutputFileName]);
info = fastqinfo(Params.ReadsFile);
len = double(info.NumberOfEntries);
Refs = fastaread(Params.RefFile);
fidOut = fopen(Params.OutputFileName,'w');
if (fidOut == -1)
    error(['Cannot create output file: ' Params.OutputFileName]);
end

% Convert all C which is not CG to T
RefsConv = ConvertNonCGtoT(Refs);


Refs.Sequence = upper(Refs.Sequence);
CGinds = findstr(Refs.Sequence,'CG');
Cinds = find(Refs.Sequence=='C');
NonCGinds = setdiff(Cinds,CGinds);
[~, CGloc] = intersect(Cinds,CGinds);
NonCGloc = setdiff(1:length(Cinds),CGloc);
Nhist = zeros(5,length(Cinds));


Nhist_CG_Reads = struct;

ReadsCount = 0;

ConvTable('A') = 1;
ConvTable('a') = 1;
ConvTable('C') = 2;
ConvTable('c') = 2;
ConvTable('G') = 3;
ConvTable('g') = 3;
ConvTable('T') = 4;
ConvTable('t') = 4;
ConvTable('-') = 5; % No Match
ConvTable('N') = 5; % No Match

fidReads = fopen(Params.ReadsFile,'r');
if (fidReads == -1)
    error(['Cannot open reads file: ' Params.ReadsFile]);
end

i=0;
line=1;
tic;
while (line ~= -1)
    i=i+1;
    if (mod(i,1000)==0)
        t = toc;
        rem = t*(len-i)/i;
        disp([name_core ': Finished ' num2str(i) '/' num2str(len) '(' num2str(round(100*i/len)) '%) Remaining: ' datestr(datenum(0,0,0,0,0,rem),'HH:MM:SS')]);

    end
    
    lineName = fgetl(fidReads);
    line = fgetl(fidReads);
    lineTmp = fgetl(fidReads);
    lineTmp = fgetl(fidReads);  
    if (lineTmp == -1)
        break;
    end    
    

    % Print read name
    fprintf(fidOut,'%s\t', lineName);

    % Ignore short reads
    if (length(line) < Params.Thresholds.MinReadLength)
        fprintf(fidOut,'Read_Length<%d\t',Params.Thresholds.MinReadLength);
        % Print original line
        fprintf(fidOut,'%s\n', line);
        continue;
    end
    foundRef = 0;

      
    [~, Al, St] = swalign(RefsConv.Sequence,line);
    % Ignore short alignments
    if (length(Al(2,:)) < Params.Thresholds.MinAlignmentLength)
        fprintf(fidOut,'Align_Length<%d\t',Params.Thresholds.MinAlignmentLength);
        % Print original line
        fprintf(fidOut,'%s\n', line);                     
        continue;
    end
    RefScore = sum(Al(2,:)=='|')/length(Al(2,:));

    % Check correct target read
    if (RefScore >= Params.Thresholds.PercentCorrectAlignment)
        ReadCheck = Al(3,find(Al(1,:)~='-'));
        CindsCheck = Cinds(find(Cinds>=St(1)))-St(1)+1;
        iShift = length(Cinds) - length(CindsCheck)+1;
        CindsCheck = CindsCheck(find(CindsCheck <= length(ReadCheck)));

        foundRef = 1;
        inds = sub2ind(size(Nhist), ConvTable(ReadCheck(CindsCheck)), iShift:length(ReadCheck(CindsCheck))+iShift-1);

        alignStr=repmat('-',1,length(RefsConv.Sequence));
        alignStr(St(1):St(1)+length(ReadCheck)-1) = ReadCheck;

        ReadsCount = ReadsCount+1;
        Nhist(inds) = Nhist(inds)+1;
        histf = alignStr(CGinds);
        histf(histf == '-') = 'X';
        if (isfield(Nhist_CG_Reads,histf) == 0)
            Nhist_CG_Reads.(histf) = 1;
        else
            Nhist_CG_Reads.(histf) = Nhist_CG_Reads.(histf)+1;
        end
    end
    if (foundRef == 0)
        fprintf(fidOut,'No_Ref_Found\t');
        % Print original line
        fprintf(fidOut,'%s\n', line);
    else   
        % Print Reference name
        fprintf(fidOut,'%s\t', Refs.Header);
        % Print original line
        fprintf(fidOut,'%s\t', line);
        % print alignment sring
        fprintf(fidOut,'%s\t', alignStr);    
        % print CG indices
        fprintf(fidOut,'CpGs:\t');
        for k=1:length(CGinds)
            fprintf(fidOut,'%s ', alignStr(CGinds(k)));    
        end
        fprintf(fidOut,'\n');    
    end
   
end

fclose(fidOut);
fclose(fidReads);

% Write summary data
Params.OutputSummary = strrep(Params.RefFile,'.fasta','.meth_summary');

fid = fopen(Params.OutputSummary,'w');
if (fid == -1)
    error(['Cannot open output summary file: ' Params.OutputSummary]);
end

fprintf(fid,'Reads file name:\t%s\n',Params.ReadsFile);
fprintf(fid,'Total Reads:\t%d\n',i);

fprintf(fid,'Sample name:\t%s\n',Refs.Header);
fprintf(fid,'# of reads aligned:\t%d\n',ReadsCount);
fprintf(fid,'CpG sites\nindex\t#A\t#C\t#G\t#T\t#-\t%%Meth\n');
for j=1:length(CGloc)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n',CGinds(j),Nhist(1,CGloc(j)),Nhist(2,CGloc(j)),Nhist(3,CGloc(j)),Nhist(4,CGloc(j)),Nhist(5,CGloc(j)),Nhist(2,CGloc(j))/(Nhist(2,CGloc(j))+Nhist(4,CGloc(j))));
end
fprintf(fid,'Non CpG sites\nindex\t#A\t#C\t#G\t#T\t#-\t%%Meth\n');
for j=1:length(NonCGloc)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n',NonCGinds(j),Nhist(1,NonCGloc(j)),Nhist(2,NonCGloc(j)),Nhist(3,NonCGloc(j)),Nhist(4,NonCGloc(j)),Nhist(5,NonCGloc(j)),Nhist(2,NonCGloc(j))/(Nhist(2,NonCGloc(j))+Nhist(4,NonCGloc(j))));
end


fclose(fid);

% Write histogram data

% CpG hist data
Params.OutputHist = strrep(Params.RefFile,'.fasta','.meth_CpG_hist');

fid = fopen(Params.OutputHist,'w');
if (fid == -1)
    error(['Cannot open output summary file: ' Params.OutputHist]);
end

for i=1:length(Refs)
    fprintf(fid,'Sample name:\t%s\n',Refs.Header);
    fprintf(fid,'Count\t');
    fprintf(fid,'# of A\t# of C\t# of G\t# of T\t# of Sites\t');
    for j=1:length(CGinds)
        fprintf(fid,'CpG Index %d\t',CGinds(j));
    end
    
    fprintf(fid,'\n');
    [a, b] = sort(struct2array(Nhist_CG_Reads),'descend');
    names = fieldnames(Nhist_CG_Reads);
    for j=1:length(a)
        name = names{b(j)};
        name(name == 'X') = '-';
        fprintf(fid,'%d\t',a(j));
        fprintf(fid,'%d\t%d\t%d\t%d\t%d\t',sum(name=='A'),sum(name=='C'),sum(name=='G'),sum(name=='T'),length(name)-sum(name=='-'));
        for k=1:length(name)
            fprintf(fid,'%c\t',name(k));
        end
        fprintf(fid,'\n');
    end
end

fclose(fid);

end
function [RefsConv] = ConvertNonCGtoT(Refs)
%% Converts CH nucleotides in sequence to TH
CGinds = findstr(Refs.Sequence,'CG');
Cinds = findstr(Refs.Sequence,'C');
RefsConv.Sequence = Refs.Sequence;
RefsConv.Sequence(setdiff(Cinds,CGinds)) = 'T';


end
function results=numToutput(cpgs,sample_name1)
%% Creates final summary table of all samples
%% Counts # of molecules with all CpGs unmethylated, all-1 CpGs unmethylated, etc.
results=zeros(1,cpgs+3);

file= [sample_name1 '.meth_CpG_hist'];

fid=fopen(file,'r');
data=textscan(fid,'%f %f %f %f %f %f %*[^\n]', 'delimiter','\t','HeaderLines', 2, 'CollectOutput',1);
fclose(fid); 

% remove irrelevant lines (not all cpgs/contains A or G)
data{1}(data{1}(:,6)<cpgs,:) = [];
data{1}(data{1}(:,2)>0,:) = [];
data{1}(data{1}(:,4)>0,:) = [];

%output2
results(1) = cpgs;

%Total aligned
results(2)=sum(data{1}(:,1));


for i = 0:cpgs
    idx=find(data{1}(:,5)==(cpgs-i));
    if ~isempty(idx)
        results(3+i) = sum(data{1}(idx,1));
    end
end

end
