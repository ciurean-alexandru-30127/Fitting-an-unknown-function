clear
load('proj_fit_39.mat');
l=length(id.X{1}); %stocam in l lungimea vectorului de identificare
l2=length(val.X{1}); %stocam in l2 lungimea vectorului de validare
MSE=[0;0]; %am creat un vector gol pentru a stoca MSE-ul in el 
for m=1:30 %m reprezinta gradul polinomului
a=[0 0]; %in a vom stoca puterile coeficientilor din aproximator
p1=0; p2=0; k=0; pozphi=1; %p1 si p2 reprezinta  puterile coeficientilor, care la final vor fi puse in a
    for p1=0:m
    		for p2=0:m
        			if p1+p2<=m
            			a=cat(1,a,[p1 p2]); %am verificat sa nu luam puteri diferite fata de gradul de care avem nevoie prin concatenare cu ajutorul functiei cat
        			end
            end
    end
p1=a(1:end,1)';
p2=a(1:end,2)'; %am pus in p1 si p2 puterile calculate ale polinomului
for i=1:l; %cu i si j am parcurs matricea de identificare
     for j=1:l 
        			for k=1:length(p1) %cu ajutorul lui k am luat pe rand puterile pentru a forma polinomul
PHY(pozphi,k)=id.X{1}(j)^p1(k)*id.X{2}(i)^p2(k); %aceasta este forma generala de calcul al elementelor matricii phy
                    end
pozphi=pozphi+1; %pozphi reprezinta numarulu de linii a matricii phy
    end
end
pozphi=pozphi-1; 
Y=reshape(id.Y',[],1); %Y fiind inital o matrice, am organizat-o intr-un vector coloana cu ajutorul functiei reshape pentru a putea calcula THETA
THETA=PHY\Y; %THETA reprezinta aproximatorul polinomial  pe care l-am aflat prin impartirea polinomiala la stanga '\'
yhatid=PHY*THETA; %yhatid l-am calculat pentru a ne ajuta la gasirea MSE
MSE(1,m)=1/pozphi*(sum((Y'-yhatid').^2)); %am aplicat formula pentru calculul MSE pentru datele de identificare
if m==6 %pentru calculul lui yhatidideal am ales m=6 pentru ca MSE pe datele de validare este minim
yhatidideal=PHY*THETA; %am calculat yhat ideal pentru gradul m=6
end

%in partea a doua, pana la partea de afisare, codul este exact acelasi doar ca este realizat pentru datele de validare
pozphi2=1;
for i=1:l2;
   		for j=1:l2 
                    for k=1:length(p1) 
PHYVAL(pozphi2,k)=val.X{1}(j)^p1(k)*val.X{2}(i)^p2(k);
        			end
        			pozphi2=pozphi2+1;
    		end
end
pozphi2=pozphi2-1;
 	YVAL=reshape(val.Y',[],1);
yhatval=PHYVAL*THETA;
MSE(2,m)=1/pozphi2*sum((YVAL'-yhatval').^2);
if m==6
    		yhatvalideal=PHYVAL*THETA;
end
end
 
plot(1:m,MSE(1,1:m),'b'),title('MSE pentru datele de identificare si pentru datele de validare')
hold on
plot(1:m,MSE(2,1:m),'r') %afisarea MSE sub forma de grafic
legend('MSE date identificare','MSE date de validare'), xlabel('Gradul polinomului: m');
ylabel('MSE');

 %reprezentarea grafica a aproximarii
figure
yhatidideal=reshape(yhatidideal,l,l)';
mesh(id.X{1},id.X{2},id.Y'), title('Datele de identificare'),xlabel('X1'),ylabel('X2'),zlabel('Yid')
figure
mesh(id.X{1},id.X{2},yhatidideal'),title('Aproximatorul polinomial ideal pentru datele de identificare'),xlabel('X1'),ylabel('X2'),zlabel('Yhatidideal')

  
figure
yhatvalideal=reshape(yhatvalideal,l2,l2)';
mesh(val.X{1},val.X{2},val.Y'),title('Datele de validare'),xlabel('X1'),ylabel('X2'),zlabel('Y')
figure
mesh(val.X{1},val.X{2},yhatvalideal'),title('Aproximatorul polinomial ideal pentru datele de validare'),xlabel('X1'),ylabel('X2'),zlabel('Yhatvalideal')

