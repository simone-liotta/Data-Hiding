%% LEGGI IMMAGINE
img = imread('gara.bmp'); %immagine di test multipla di 512x512
%% PARAMETRI DI PROGETTO
alfa = 5; %alfa nella marchiatura ADDITIVA di Cox
peso = 3; %maggiorazione di alfa per aree a varianza alta
varTresh = 0.8e+03; %soglia di varianza alta
w = randi([0,1],1,1024); %marchio 

%% ANALIZZO L'IMMAGINE A BLOCCHI DI 512x512
fun = @(block_struct) GR4embedding(block_struct.data,w,alfa,peso,varTresh);
I2 = blockproc(img,[512 512],fun,'PadPartialBlocks',true);

% Salva immagine marchiata su disco
imwrite(uint8(I2),'garaMarchio.bmp')
%% SE L'IMMAGINE NON FOSSE STATA MULTIPLO DI 512x512
%     LG = size(Ig) ~= size(I2);
%     if ismember(1,LG) == 1
%         
%         r = size(Ig,1);
%         c = size(Ig,2);
%         
%         I2 = I2(1:r,1:c);
%         
%     end
%% CALCOLO IL WPSNR DELL'IMMAGINE MARCHIATA RISPETTO ALL'ORIGINALE
Id = double(img); 
Iwmax = I2./max(max(I2));
Imax = Id./max(max(Id));
wpsnr_out = WPSNR(abs(Imax),abs(Iwmax))
        
figure,imshow(uint8(I2))


