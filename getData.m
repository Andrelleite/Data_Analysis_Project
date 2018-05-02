
%------INDICE DE FUNÇÕES UTILIZADAS------

%Andre Correia 
%Andre Leite

% isnan(vect) -> verifica se tem valores NAN no vetor;
% mean(vect) -> calcula a média
% std(vect) -> calcula o desvio padrão
% detrend(x) -> remove o melhor ajuste no vetor x e faz o seu retorno(tendencia)
% detrend (x , 'constant') -> remove o valor médio do vetor x
% detrend (x , 'linear') -> remove a tendencia linear do vetor x
% polyfit(x,y,n) -> retorna os coeficientes para uma polinomial de grau n que melhor se adequa à informação de y.
% polyval(p,x) -> retorna o valor de uma polinomial de grau n testada em x
% repmat(A,r1,...,rN) -> specifica uma lista de escalares, r1,..,rN, que descrever como as copias de A são organizados em cada dimensão.
% dummyvar(temp) -> retorna a matriz D contendo zeros e uns, cujas colunas são varáveis "dummy" para temp
% x = A\B -> A * x = B
% adftest(vect) -> retorna um valor com a decisão de rejeição
% autocorr(y1) -> indicador para a ordem de nc(modelo MA)
% parcorr(y1) -> indicador para a ordem de na(modelo AR)

%----------------------------------------

function [dados,tempo] = getData()
    
    disp('PROJETO - Análise e Transformação de Dados - 2018');
    fprintf('\nNumero da sua PL (pratica laboratorial)\n');
    PL = input( 'PL: ' , 's');
    fich = "dataset_ATD_PL"+PL+".csv";
    disp(fich);
    
    allData = load(fich); %load dos dados do dataset (nao pode ter cabeçalho nas tabelas)
    time = allData(1:end,1); %obtenção de dados na primeira coluna
    data = allData(1:end,2); %obtenção de dados na segunda coluna
    dados = data; 
    tempo = time;
    treatData(dados); %passagem de dados para tratamento em função

end
 
function treatData(x1)
           
    N=length(x1); %comprimento da serie temporais
    t=(0:N-1)'; %escala temporal
    x1r=x1; %atribuição a variavel temporária

    %________________________________1_____________________________________
    
    if any(isnan(x1)) %pesquisa de valores NAN
        ind=find(isnan(x1)); %indices de NAN em vetor
        for k=1:length(ind) 
            tt=t(ind(k)-2:ind(k)-1);
            xx=x1r(ind(k)-2:ind(k)-1);
            x1r(ind(k))=interp1(tt,xx,t(ind(k)),'pchip','extrap'); %interpolação com método pchip para substituição de valores NAN
        end
    end
    
    x1n = x1r;
    disp(x1n); %visualização de vetor com NANs preenchidos
    
    media = mean(x1r); %calculo da media 
    padrao = std(x1r); %calculo do desvio padrao
    max = media + (1*padrao); %valor para secção mensal (+desvio)
    min = media - (1*padrao); %valor para secção mensal (-desvio)

    for i=1:length(x1r) %check de outliers
        if x1r(i)>max %para valores acima da média + desvio padrão
            x1r(i)= media+(2.5*padrao); %resolução de outliers (LOW)
        elseif x1r<min %para valores abaixo da média - desvio padrão
            x1r(i)= media-(2.5*padrao); %resolução de outliers (HIGH)
        end
    end
    
    %Sumarização gráfica e comparação de dados adquiridos
    figure(1)
    subplot(311) 
    plot(t,x1,'+-');
    legend('serie temporal Original','Location','northwest')
    xlabel('t[h]');
    title('series temporais');
    subplot(312)
    plot(t,x1n,'+-');
    legend('serie temporal sem NANs e com outliers','Location','northwest')
    xlabel('t[h]');
    subplot(313)
    plot(t,x1r,'+-');
    legend('serie temporal sem NANs e sem outliers','Location','northwest')
    xlabel('t[h]');
       
    %_____________________________2________________________________________
    
    x1rz = detrend(x1r,'constant'); %aproximação polinomial de ordem 0
    trz = x1r - x1rz; %tendencia de grau 0
    x1ru = detrend(x1r,'linear'); %aproximação polinomial de ordem 1
    tru = x1r - x1ru; %tendencia de grau 1
    
    %Sumarização gráfica e comparação de dados
    figure(2)
    subplot(621)
    plot(t,x1r,'+-r',t,trz,'+-b');
    legend('Serie regularizada','Tendencia de ordem 0','Location','northwest');
    subplot(622)
    plot(t,x1rz,'+-r');
    legend('Serie sem tendencia de ordem 0','Location','northwest');
    subplot(623)
    plot(t,x1r,'+-r',t,tru,'+-b');
    legend('Serie regularizada','Tendencia de ordem 1','Location','northwest');
    subplot(624)
    plot(t,x1ru,'+-r');
    legend('Serie sem tendencia de ordem 1','Location','northwest');
    subplot(625)
    plot(t,trz,'+-g');
    legend('Tendencia de ordem 0','Location','northwest');
    subplot(626)
    plot(t,tru,'+-r');
    legend('Tendencia de ordem 1','Location','northwest');
    
    pl1 = polyfit(t,x1r,2); %coeficientes dos polinomios
    trd = polyval(pl1,t); %tendencia de ordem 3
    x1rd = x1r - trd; %calculo de serie reguralizada sem tendencia
    
    %Sumarização Gráfica e comparação de dados
    figure(3)
    subplot(331)
    plot(t,x1r,'+-r',t,trd,'+-g');
    legend('Serie regularizada','Tendencia de ordem 3','Location','northwest');
    xlabel('t[h]');
    subplot(332)
    plot(t,x1rd,'+-r');
    legend('Serie sem tendencia de ordem 3','Location','northwest');
    xlabel('t[h]');
    subplot(333)
    plot(t,trd,'+-g');
    legend('Tendencia de ordem 3','Location','northwest');
    xlabel('t[h]');
    
    %Estimativa Sazonal
    temp = repmat([1:7]',52,1); %sazonalidade semanal
    sz = dummyvar(temp); %gera matriz dummy
    x1r = x1r(1:end-1,:);
    x1rd = x1rd(1:end-1,:);
    bs = sz\x1rd;
    st = sz*bs; %componente sazonal
    x1s = x1rd - st;
    t = t(1:end-1,:);
    
    %Sumarização Gráfica e comparação de dados
    figure(4)
    subplot(421)
    plot(t,x1rd,'+-r',t,st,'+-g');
    legend('Serie sem tendencia de ordem 4','Sazonalidade','Location','northwest');
    xlabel('t[h]');
    subplot(422)
    plot(t,x1s,'+-r');
    legend('Serie sem Sazonalidade','Location','northwest');
    xlabel('t[h]');
    
    irreg = x1r-pl1-st; %componente irregular 
    snirreg = x1r-irreg; %serie temporal sem componente irregular
    
    %Sumarização Gráfica e comparação de dados    
    figure(5)
    subplot(521)
    plot(t,snirreg,'+-r');
    legend('Serie sem componente irregular','Location','northwest');
    subplot(522)
    plot(t,irreg,'+-r');
    legend('Serie com componente irregular','Location','northwest');    
    
    %_______________________________3______________________________________
    
    tt=(0:2*N-1)'; %escala temporal para previsão

    
    adftest(x1r); %teste de estacionaridade da série regularizada
    adftest(st); %teste de estacionaridade da componente sazonal
    
    t1=(0:6)';
    y1=st(1:7); %componente sazonal
    disp(y1);
    
    figure(6)
    subplot(621)
    autocorr(y1); %indicador para a ordem de nc(modelo MA)
    subplot(622)
    parcorr(y1); %indicador para a ordem de na(modelo AR)
        
    id_y1 = iddata(y1, [], 1, 'TimeUnit', 'days');
    
    opt1_AR = arOptions('Approach', 'ls');
    na1_AR = 3; %histórico da variável
    model1_AR = ar(id_y1,na1_AR, opt1_AR); %modelo AR
    pcoef1_AR = polydata(model1_AR); %parâmetros do modelo AR
    
    y1_AR = y1(1:na1_AR);
    
    for k=na1_AR+1:7
        y1_AR(k)=sum(-pcoef1_AR(2:end)'.*flip(y1_AR(k-na1_AR:k-1)));
    end
    
    y1_AR2=repmat(y1_AR,52,1);
    
    y1_ARf=forecast(model1_AR,y1(1:na1_AR),7-na1_AR);
    y1_ARf2=repmat([y1(1:na1_AR); y1_ARf],52,1);
    
    figure(7)
    plot(t,st,'-+',t,y1_AR2,'-o',t,y1_ARf2,'-*');
    xlabel('t [h]');
    title('Componente sazonal 1 (-+) e estimação com modelo AR');
    figure(8)
    plot(t,x1r,'-+',t,y1_AR2, t, trd(1:end-1) ,'-o');
    xlabel('t [h]');
    title('Série 1(-+) e estimação com o modelo AR(-o)');
        
    E1_AR = sum((x1r-(y1_AR2+trd(1:end-1))).^2);
    
    tr1_2_2=polyval(pl1,tt); %tendencia

%     figure(9)
%     plot(t,x1rd,'-+',tt,repmat(y1_AR2,52,1)+tr1_2_2,'-o');
%     xlabel('t [h]');
%     title('Série 1 (-+) e Previsão com o modelo AR (-o)');
    
    % Estimação de um modelo arma
    opt1_ARMAX = armaxOptions('SearchMethod', 'auto');
    na1_ARMA=1;
    nc1_ARMA=1;
    model1_ARMA = armax(id_y1,[na1_ARMA nc1_ARMA], opt1_ARMAX);
    [pa1_ARMA,pb1_ARMA,pc1_ARMA] = polydata(model1_ARMA);
    
    e = randn(7,1); %ruído branco
    y1_ARMA = y1(1:na1_ARMA);
    
    for k=na1_ARMA+1:7
        y1_ARMA(k)=sum(-pa1_ARMA(2:end)'.*flip(y1_ARMA(k-na1_ARMA:k-1)))+sum(pc1_ARMA'.*flip(e(k-nc1_ARMA:k)));
    end
    y1_ARMA2=repmat(y1_ARMA,52,1);

    %Simulação do modelo arma com forecast
    y1_ARMAf=forecast(model1_ARMA, y1(1:na1_ARMA),7-na1_ARMA);
    y1_ARMAf2=repmat([y1(1:na1_ARMA); y1_ARMAf],52,1);
    y1_ARMA2 = y1_ARMA2(:);
    figure(10) %compara a componente sazonal com a sua estimação
    plot(t,st,'-+',t,y1_ARMA2,'-o',t,y1_ARMAf2,'-*');
    xlabel('t [h]');
    title('Componente sazonal 1(-+) e estimação com o modelo ARMA');

    figure(11) %compara a série com o modelo ARMA + tendência
    plot(t,x1r,'-+',t,y1_ARMA2+trd(1:end-1),'-o');
    xlabel('t [h]');
    title('Série 1 (-+) e estimação com o modelo ARMA(-o)')

    %Métrica para análise
    E1_ARMA = sum((x1r-y1_ARMA2(1:N-1)).^2);

%     figure(12) %Faz a previsão para 2N
%     plot(t,x1r,'-+',tt,repmat(y1_ARMA2,52,1)+tr1_2_2, '-o');
%     xlabel('t [h]');
%     title('Série 1(-+) e previsão com o modelo ARMA(-o)')

    %Estimação de um modelo ARIMA
    p1_ARIMA = 3;
    q1_ARIMA = 3;
    Mdl = arima(p1_ARIMA,0,q1_ARIMA);
    EstMd1 = estimate(Mdl,x1r(1:N-1),"Y0",x1r(1:p1_ARIMA+1));
    
    %Simulação do modelo ARIMA
    y1_ARIMA = simulate(EstMd1,N-1);

    figure(13) %compara a serie com a sua estimaçao
    plot(t,x1r,'-+',t,y1_ARIMA,'-o');
    xlabel('t [h]');
    title('Série 1(-+) e estimação com o modelo ARIMA(-o)')

    %métrica para análise
    E1_ARIMA = sum((x1r-y1_ARIMA(1:N-1)).^2);
    
    %Simulação do modelo ARIMA para 2N
    y1_ARIMA2 = simulate(EstMd1,2*N);

    figure(14) %faz a previsão para 2N
    plot(t,x1r,'-+',tt,y1_ARIMA2,'-o');
    xlabel('t [h]');
    title('Série 1(-+) e estimação com o modelo ARIMA(-o)')
    
end


