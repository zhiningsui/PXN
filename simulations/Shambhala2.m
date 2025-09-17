fileID = fopen('args.txt','r'); % arg.txt contains the number of samples to be normalized, the number of samples in auxiliary matrix, and k. 
sizeA = [3 1];
A = fscanf(fileID,'%f',sizeA); % "A" contains the content of arg.txt
fclose(fileID); 
    
NH = A(1); % the number of samples to be normalized
NP = A(2); % the number of samples in the pooled matrix P
k = A(3);
 
inData=readExpressionData("P_prim.txt", "not_log2"); 
Exp = inData.Samples;
SYMBOL = inData.GeneList;
SN = inData.SamplesName;
NS = length(SN);

% Harmonize samples one-by-one
for ( i = 1:NH ) 
    message = sprintf('Harmonizing sample %d out of %d',i,NH);
    disp(message);
    
    % i0: the column indices for the sample i and all samples in the
    % auxiliary matrix.
    i0 = i;
    for ( jjj = 1:NP )  % loop over samples in the auxiliary matrix          
        i0 = [i0 (jjj+NH)];
    end   
        
    EXP = Exp(:,i0); % the expression data for sample i and auxiliary samples
    EXP = quantilenorm(EXP); % quantile-normalize EXP
        
    dataN = CuBlock(real(EXP),[],k); % Run CuBlock to normalize EXP  

    log2e = log2(exp(1));
    DataN = dataN/log2e;

    EXPN = exp(DataN)-1;
    vecN = EXPN(:,1);
        
    if ( i == 1 ) 
        OUT = vecN;
    else    
        OUT = [OUT vecN];
    end

end

outFileName = append("Cu_bis.txt"); 
outFile=fopen(outFileName,'w');
nG=size(OUT,1);
nS=size(OUT,2);
for i=1:nG
    fprintf(outFile,'%s',SYMBOL{i,1});        
    for j=1:nS
        fprintf(outFile,' %f',OUT(i,j));        
    end
    fprintf(outFile, '\n');
end
fclose(outFile);
    

