function [Y,HL14,HL24] = LIOTTAembedding(img,w,alfa,peso,varTresh)
%% PASSO NELLO SPAZIO DI COLORI 'YST' (se considero l'analisi di volti)
%Yst_coeff = [1 0 0; 0 0.7072 0.6345; 0 -0.6345 0.7072];

% r = double(imgRGB(:,:,1));
% g = double(imgRGB(:,:,2));
% b = double(imgRGB(:,:,3));


% Estraggo le componenti YUV da definizione (vedi wikipedia)
% Wr = 0.299;
% Wb = 0.114;
% Wg = 0.587;
% Umax = 0.436;
% Vmax = 0.615;
%
% Y = Wr*r+Wg*g+Wb*b;
%
% U = Umax*((b-Y)/(1-Wb));
%
% V = Vmax*((r-Y)/(1-Wr));


% Componenti dominio YST (Neri - Giunta - Benedetto)
% S = 0 + 0.7072*U + 0.6345*V; %componente Skin (S)
% %imshow(S)
% T = 0 - 0.6345*U + 0.7072*V; %componente ortogonale a 'skin' (T)
%imshow(T)
%% ANALISI AUTOMATICA DELLA VARIANZA DELL'IMMAGINE

n = 64*ones(1,8);
Ymat = mat2cell(img,n,n); %creo una memoria contenente 8 blocchi 64x64
Yvett = reshape(Ymat',1,[]);

for k = 1:64
    yvInt = Yvett{1,k};
    v(1,k) = var(double(yvInt(:))); %vettore con varianze memorizzate
end


vAlta = find(v > varTresh); %seleziona solo i blocchi con varianza maggiore di soglia

%% APPLICO L'ALGORITMO 'DIGITAL IMAGE WATERMARKING BASED ON DWT-DCT'
% DWT A 4 LIVELLI

% STEP 1 (sottobande 256x256)
[LL1,HL1,LH1,HH1] = dwt2(img,'haar'); % dwt su intera immagine
% imshow(uint8(LL1))
% imshow(uint8(HL1))
% imshow(uint8(LH1))
% imshow(uint8(HH1))

% STEP 2 (sottobande 128x128)
[LL12,HL12,LH12,HH12] = dwt2(HL1,'haar'); % dwt su componente HL1
[LL22,HL22,LH22,HH22] = dwt2(LH1,'haar'); % dwt su componente LH1

% STEP 3 (sottobande 64x64)
[LL13,HL13,LH13,HH13] = dwt2(HL12,'haar'); % dwt su componente HL12
[LL23,HL23,LH23,HH23] = dwt2(LH22,'haar'); % dwt su componente LH22

% STEP 4 (sottobande 32x32)
[LL14,HL14,LH14,HH14] = dwt2(HL13,'haar'); % dwt su componente HL13
[LL24,HL24,LH24,HH24] = dwt2(LH23,'haar'); % dwt su componente LH23

%% MARCHIO LA SOTTOBANDA HL14
% divido le componenti HL14 in blocchi 4x4
N = 4*ones(1,8);
blockHL14 = mat2cell(HL14,N,N); %creo una memoria contenente 64 blocchi 4x4
vettHL24 = reshape(blockHL14',1,[]); %per semplicità riscrivo la matrice come un vettore riga

%riordino il vettore del marchio nella forma delle componenti HL14 e LH24
wMat = reshape(w,32,32).'; %riscrivo il vettore del rumore come una matrice
blockW = mat2cell(wMat,N,N); %divido la matrice del rumore in blocchi 4x4
vettW = reshape(blockW',1,[]); %riscrivo di nuovo la matrice con i blocchi in riga per semplicità

%Performo la DCT su tutti e 64 i blocchi 4x4 e amrchio con regola additiva
%di Cox
for i=1:64 
    imgDCT = dct2(vettHL24{1,i}); %fai la dct del blocco i-esimo di sottobanda di immagine
    wDCT = dct2(vettW{1,i});%fai la dct del blocco i-esimo di rumore
    
    trovato = ismember(i,vAlta);% seleziona i blocchi con varianza alta
    
    if trovato == 1
        % nei blocchia  varianza alta marchi con un alfa maggiore
        vettHL24{1,i} = idct2(imgDCT(:,:)+ ((alfa+peso)*wDCT(:,:)) );
    else
        %COX'S ADDITIVE ALGORITHM EMBEDDING
        vettHL24{1,i} = idct2(imgDCT(:,:)+(alfa*wDCT(:,:)) );
    end    
end

HL14 = cell2mat(reshape(vettHL24,8,8).'); %restituzione della sottobanda HL14 marchiata

%% MARCHIO LA SOTTOBANDA LH24
% divido le componenti LH14 in blocchi 4x4
blockLH24 = mat2cell(LH24,N,N); %creo una memoria contenente 64 blocchi 4x4
vettLH24 = reshape(blockLH24',1,[]); %per semplicità riscrivo la matrice come un vettore riga

%Performo la DCT su tutti e 64 i blocchi 4x4
for j=1:64
    %iteraBlocco = vettHL14{1,i};
    imgDCT = dct2(vettLH24{1,j});
    wDCT = dct2(vettW{1,j});
    
    trovato = ismember(j,vAlta);% seleziona i blocchi con varianza alta
    
    if trovato == 1
        % nei blocchi a  varianza alta marchi con un alfa maggiore
        vettHL24{1,j} = idct2(imgDCT(:,:)+ ((alfa+peso)*wDCT(:,:)) );
    else
        %COX'S ADDITIVE ALGORITHM EMBEDDING
        vettHL24{1,j} = idct2(imgDCT(:,:)+ (alfa*wDCT(:,:)) );   
    end
end

LH24 = cell2mat(reshape(vettLH24,8,8).'); %restituzione della sottobanda HL14 marchiata

%% RICOSTRUISCO L'IMMAGINE MARCHIATA

%ricostruisco la sottobanda 64x64
HL13 = idwt2(LL14,HL14,LH14,HH14,'haar');
LH23 = idwt2(LL24,HL24,LH24,HH24,'haar');

%ricostruisco la sottobanda 128x128
HL12 = idwt2(LL13,HL13,LH13,HH13,'haar');
LH22 = idwt2(LL23,HL23,LH23,HH23,'haar');

%ricostruisco la sottobanda 256x256
HL1 = idwt2(LL12,HL12,LH12,HH12,'haar');
LH1 = idwt2(LL22,HL22,LH22,HH22,'haar');

%ricostruisco l'intera componente T MARCHIATA 512x512
Tw = idwt2(LL1,HL1,LH1,HH1,'haar');
Y = Tw; %Se marchio sulla componente della luminanza Y (per i NON volti)
%% RITORNO AL DOMINIO RGB da YST (Neri - Giunta - Benedetto)

% U = 0.7834*S - (0.7029*((Tw))); %ritorno a componente U da YST marchiato
% V = 0.7029*S + (0.7834*((Tw))); %ritorno a componente V da YST marchiato

%torno in RGB da YUV
% r = Y + V*((1-Wr)/Vmax);
% g = Y - (U*((Wb*(1-Wb))/Umax*Wg))- (V*((Wr*(1-Wr))/Vmax*Wg));
% b = Y + U*((1-Wb)/Umax);

% Iw = cat(3,r,g,b);
%ricostruisco l'immagine originale marchiata dalle componenti RGB
% Iw = zeros(size(imgRGB));
% Iw(:,:,1) = r;
% Iw(:,:,2) = g;
% Iw(:,:,3) = b;
%figure; imshow(uint8(Iw))
end