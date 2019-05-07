% PASSA IN INPUT (COME STRINGA) IL 'nome.estensione'
% DELL'IMMAGINE ORIGINALE E QUELLO DELL'IMMAGINE MARCHIATA

function [MarchioPresente] = LIOTTAdetection(nomeImOrig, nomeImMarchiata,w)
%% LEGGO LE IMMAGINI E SETTO I PARAMETRI DI PROGETTO
imgO=imread(nomeImOrig);
imgM=imread(nomeImMarchiata);

alfa = 5;
peso = 3;
varTresh = 0.8e+03;

%% DIVIDO IN 4 ZONE L'IMMAGINE MARCHIATA
Nb = 512*ones(1,2);
matZonaM = mat2cell(imgM,Nb,Nb); %creo una memoria contenente 4 blocchi 512x512
ZonaM = reshape(matZonaM',1,[]); %per semplicità riscrivo la matrice come un vettore riga

%% DIVIDO IN 4 ZONE L'IMMAGINE ORIGINALE
matZonaO = mat2cell(imgO,Nb,Nb); %creo una memoria contenente 4 blocchi 512x512
ZonaO = reshape(matZonaO',1,[]); %per semplicità riscrivo la matrice come un vettore riga

%% ELABORO CIASCUNA DELLE 4 AREE 512x512
for zona = 1 : 4
    %% CALCOLO DWT2 PER L'IMMAGINE ORIGINALE
    [LL1o,HL1o,LH1o,HH1o] = dwt2(ZonaO{1,zona},'haar'); % dwt su intera immagine
    
    % STEP 2 (sottobande 128x128)
    [LL12o,HL12o,LH12o,HH12o] = dwt2(HL1o,'haar'); % dwt su componente HL1
    [LL22o,HL22o,LH22o,HH22o] = dwt2(LH1o,'haar'); % dwt su componente LH1
    
    % STEP 3 (sottobande 64x64)
    [LL13o,HL13o,LH13o,HH13o] = dwt2(HL12o,'haar'); % dwt su componente HL12
    [LL23o,HL23o,LH23o,HH23o] = dwt2(LH22o,'haar'); % dwt su componente LH22
    
    % STEP 4 (sottobande 32x32)
    [LL14o,HL14o,LH14o,HH14o] = dwt2(HL13o,'haar'); % dwt su componente HL13
    [LL24o,HL24o,LH24o,HH24o] = dwt2(LH23o,'haar'); % dwt su componente LH23
    
    %% CALCOLO DWT2 PER L'IMMAGINE MARCHIATA
    [LL1,HL1,LH1,HH1] = dwt2(ZonaM{1,zona},'haar'); % dwt su intera immagine
    
    % STEP 2 (sottobande 128x128)
    [LL12,HL12,LH12,HH12] = dwt2(HL1,'haar'); % dwt su componente HL1
    [LL22,HL22,LH22,HH22] = dwt2(LH1,'haar'); % dwt su componente LH1
    
    % STEP 3 (sottobande 64x64)
    [LL13,HL13,LH13,HH13] = dwt2(HL12,'haar'); % dwt su componente HL12
    [LL23,HL23,LH23,HH23] = dwt2(LH22,'haar'); % dwt su componente LH22
    
    % STEP 4 (sottobande 32x32)
    [LL14,HL14,LH14,HH14] = dwt2(HL13,'haar'); % dwt su componente HL13
    [LL24,HL24,LH24,HH24] = dwt2(LH23,'haar'); % dwt su componente LH23
    
    %% ANALISI AUTOMATICA DELLA VARIANZA DELL'IMMAGINE ORIGINALE
    n = 64*ones(1,8);
    Ymat = mat2cell(ZonaO{1,zona},n,n); %creo una memoria contenente 8 blocchi 64x64
    Yvett = reshape(Ymat',1,[]);
    
    for k = 1:64
        yvInt = Yvett{1,k};
        v(1,k) = var(double(yvInt(:))); %vettore con varianze memorizzate
    end
    
    vAlta = find(v > varTresh); %seleziona solo i blocchi con varianza maggiore di soglia
    
    %% ESTRAZIONE DEL MARCHIO NELLA SOTTOBANDA HL14
    % FACCIO DETECTION PER LA SOTTOBANDA HL14
    
    % CREO I BLOCCHI MARCHIATI 4X4
    N = 4*ones(1,8);
    blockHL14m = mat2cell(HL14,N,N); %creo una memoria contenente 64 blocchi 4x4
    vettHL24m = reshape(blockHL14m',1,[]); %per semplicità riscrivo la matrice come un vettore riga
    
    % CREO I BLOCCHI ORIGINALI 4X4
    N = 4*ones(1,8);
    blockHL14o = mat2cell(HL14o,N,N); %creo una memoria contenente 64 blocchi 4x4
    vettHL24o = reshape(blockHL14o',1,[]); %per semplicità riscrivo la matrice come un vettore riga
    
    
    wMat = reshape(w,32,32).'; %riscrivo il vettore del rumore come una matrice
    blockW = mat2cell(wMat,N,N); %divido la matrice del rumore in blocchi 4x4
    vettW = reshape(blockW',1,[]); %riscrivo di nuovo la matrice con i blocchi in riga per semplicità
    
    % Performo la DCT su tutti e 64 i blocchi 4x4
    for i=1:64
        
        imgMDCT = dct2(vettHL24m{1,i}); %fai la dct del blocco i-esimo di sottobanda di immagine
        imgODCT = dct2(vettHL24o{1,i}); %fai la dct del blocco i-esimo di sottobanda di immagine
        wDCT = dct2(vettW{1,i});%fai la dct del blocco i-esimo di rumore
        
        trovato = ismember(i,vAlta);% seleziona i blocchi con varianza alta
        
        if trovato == 1
            % nei blocchi a  varianza alta marchi con un alfa maggiore
            WTildeHL{1,i} = (idct2((imgMDCT(:,:)-imgODCT(:,:))/(alfa+peso)));
        else
            %COX'S ADDITIVE ALGORITHM EMBEDDING
            WTildeHL{1,i} = (idct2((imgMDCT(:,:)-imgODCT(:,:))/alfa));
        end
    end
    
    maxExtHL = max((reshape((cell2mat(reshape(WTildeHL,8,8).'))',1,1024)));
    wEstrattoHL = uint8((reshape((cell2mat(reshape(WTildeHL,8,8).'))',1,1024)))/maxExtHL;
    WtotHL{1,zona} = wEstrattoHL; %memoria contenente 4 marchi HL estratti dall'immagine    
    
    %% ESTRAZIONE DEL MARCHIO NELLA SOTTOBANDA LH24
    % FACCIO DETECTION PER LA SOTTOBANDA LH24
    
    % CREO I BLOCCHI MARCHIATI 4X4
    N = 4*ones(1,8);
    blockLH24m = mat2cell(LH24,N,N); %creo una memoria contenente 64 blocchi 4x4
    vettLH24m = reshape(blockLH24m',1,[]); %per semplicità riscrivo la matrice come un vettore riga
    
    % CREO I BLOCCHI ORIGINALI 4X4
    N = 4*ones(1,8);
    blockLH24o = mat2cell(LH24o,N,N); %creo una memoria contenente 64 blocchi 4x4
    vettLH24o = reshape(blockLH24o',1,[]); %per semplicità riscrivo la matrice come un vettore riga
    
    % Performo la DCT su tutti e 64 i blocchi 4x4
    for j=1:64
        imgMDCT = dct2(vettLH24m{1,j}); %fai la dct del blocco i-esimo di sottobanda di immagine
        imgODCT = dct2(vettLH24o{1,j}); %fai la dct del blocco i-esimo di sottobanda di immagine
        wDCT = dct2(vettW{1,j});%fai la dct del blocco i-esimo di rumore
        
        trovato = ismember(j,vAlta);% seleziona i blocchi con varianza alta
        
        if trovato == 1
            % nei blocchi a  varianza alta marchi con un alfa maggiore
            WTildeLH{1,j} = (idct2((imgMDCT(:,:)-imgODCT(:,:))/(alfa+peso)));
        else
            %COX'S ADDITIVE ALGORITHM EMBEDDING
            WTildeLH{1,j} = (idct2((imgMDCT(:,:)-imgODCT(:,:))/alfa));
        end
    end
    
    maxExtLH = max((reshape((cell2mat(reshape(WTildeLH,8,8).'))',1,1024)));
    wEstrattoLH = uint8((reshape((cell2mat(reshape(WTildeLH,8,8).'))',1,1024)))/maxExtLH;
    WtotLH{1,zona} = wEstrattoLH; %memoria contenente 4 marchi LH estratti dall'immagine
end

Wtot = cat(2,WtotHL,WtotLH);%memoria contenente gli 8 marchi estratti dall'immagine

%vado a selezionare tra gli 8 estratti il marchio più correlato a quello
%originale
maxCor=0;
indice=1;
for ind=1:8
    rho=corr2(Wtot{1,ind},w);
    if rho >maxCor
        maxCor = rho;
        indice = ind;
    end
end

%MARCHIO ESTRATTO
Westratto = Wtot{1,indice};

%CALCOLO DELLA SOGLIA
Max=0;

for u=1:1000
    Wi = randi([0,1],1,length(Westratto));
    Correlazione = corr2(Wi,w); %%%(tra il marchio rand e quello orig o estra)%%%
    if Correlazione>Max
        Max = Correlazione;
    end
end

%la soglia s è definita come:
%l'aumento del 20% del massimo della correlazione prima calcolata
s = Max+Max/5;

%operazione di confronto con la soglia, restituisce 1 per valori sopra
%soglia (presenza del marchio e 0 per i valori sottosoglia, assenza del
%marchio.
if maxCor > s
    MarchioPresente = 1;
else
    MarchioPresente = 0;
end
end
