% ce programme créé la matrice qui contient les valeurs du total de pixels
% représentant un front pour chaque mois pour 25 ans et créer le graphique
% de l'indice frontal F1 pour chacune des régions (Labrador, Bay, Groenland
% et Banks 

summary_Groenland=zeros(12,25);
summary_Bay=zeros(12,25);
summary_Banks=zeros(12,25);
summary_Labrador=zeros(12,25);

% months=['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];
% years=[1986;1987;1988;1989;1990;1991;1992;1993;1994;1995;1996;1997;1998;1999;2000;2001;2002;2003;2004;2005;2006;2007;2008;2009;2010];

source = 'H:\matlab_atlantic\subset_output\';

Dir_source = dir(source); 
dim_source = size(Dir_source); 

for i=3:dim_source(1) 
    for y=1986:2010
        year=num2str(y);
        for m=1:12
            month=m;
            if (month<10)        
                strmonth= ['0' num2str(month)];
            else
                strmonth= num2str(month);
            end
    %                                        Subset_1986_01.mat
            if (strcmp(Dir_source(i).name(8:11),year)&& strcmp(Dir_source(i).name(13:14),strmonth)) % strcmp : string comparaison
    %          s'assure que - les caractères 8 à 11 du nom de la source correspond à l'année traitée 
    %                       - les caractères 13 et 14 du nom de la source correspond au mois traitée 
    %                       - l'extension est .mat
                data=load([source Dir_source(i).name]);                       
                total_Groenland=data.total_Groenland;
                total_Bay=data.total_Bay;
                total_Banks=data.total_Banks;
                total_Labrador=data.total_Labrador;
                
                j=m;
                k=y-1985;
                summary_Groenland(j,k)=total_Groenland;
                summary_Bay(j,k)=total_Bay;
                summary_Banks(j,k)=total_Banks;
                summary_Labrador(j,k)=total_Labrador;
            end
        end
    end    
end

months=[1;2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12];
years=[1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010];

dest ='H:\matlab_atlantic\subset_output\';  
%                         
save([dest 'Groenland_fronts_1986_2010_summary.mat'],'summary_Groenland','months','years');
save([dest 'Bay_fronts_1986_2010_summary.mat'],'summary_Bay','months','years');
save([dest 'Banks_fronts_1986_2010_summary.mat'],'summary_Banks','months','years');
save([dest 'Labrador_fronts_1986_2010_summary.mat'],'summary_Labrador','months','years');


% LES FIGURES DEVRONT ÊTRE AMÉLIORÉES POUR BIEN REPRÉSENTER LES DONNÉES

 h=figure(1);
 surf(years,months,summary_Labrador)
%  axis tight
%  daspect([12 1 10])
 view (-50,30)
 colormap jet;
            
 title('Temporal variability of the frontal index F1 for the Labrador region, January 1986 - December 2010'  ) 
 
 saveas(h,[dest 'Frontal_index_Labrador_1986-2010 '],'png') 
 saveas(h,[dest 'Frontal_index_Labrador_1986-2010 '],'fig')
 delete(h)
 
 h=figure(1);
 surf(years,months,summary_Groenland)
%  axis tight
%  daspect([12 1 10])
 view (-50,30)
 colormap jet;
            
 title('Temporal variability of the frontal index F1 for the Groenland region, January 1986 - December 2010'  ) 
 
 saveas(h,[dest 'Frontal_index_Groenland_1986-2010 '],'png') 
 saveas(h,[dest 'Frontal_index_Groenland_1986-2010 '],'fig')  
 delete(h)
 
 h=figure(1);
 surf(years,months,summary_Bay)
%  axis tight
%  daspect([12 1 10])
 view (-50,30)
 colormap jet;
            
 title('Temporal variability of the frontal index F1 for the Bay region, January 1986 - December 2010'  ) 
 saveas(h,[dest 'Frontal_index_Bay_1986-2010 '],'png') 
 saveas(h,[dest 'Frontal_index_Bay_1986-2010 '],'fig') 
 
 delete(h)
 
 h=figure(1);
 surf(years,months,summary_Banks)
%  axis tight
%  daspect([12 1 10])
 view (-50,30)
 colormap jet;
            
 title('Temporal variability of the frontal index F1 for the Banks region, January 1986 - December 2010'  ) 
 saveas(h,[dest 'Frontal_index_Banks_1986-2010 '],'png') 
 saveas(h,[dest 'Frontal_index_Banks_1986-2010 '],'fig') 
 delete(h)