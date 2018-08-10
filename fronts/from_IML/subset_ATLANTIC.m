

source = 'H:\matlab_atlantic\moyenne_INPUT\';

Dir_source = dir(source); 
dim_source = size(Dir_source); 

extension  = 'mat';
     
for i=3:dim_source(1) 
    for y=1986:2010;
        year=num2str(y);
            for m=1:12    % A MODIFIER selon les mois à traiter
               month=m;
                if (month<10)        
                    strmonth= ['0' num2str(month)];
                else
                    strmonth= num2str(month);
                end
%                                           Prob_Mensuel_2010_09.mat
                if (strcmp(Dir_source(i).name(14:17),year)&& strcmp(Dir_source(i).name(19:20),strmonth)&& strcmp(Dir_source(i).name(end-2:end),extension)) % strcmp : string comparaison
        %          s'assure que - les caractères 14 à 17 du nom de la source correspond à l'année traitée 
        %                       - les caractères 19 et 20 du nom de la source correspond au mois traitée 
        %                       - l'extension est .mat
                    data=load([source Dir_source(i).name]);                       
                    original=data.count_edge;   

                    a=1:672;
                    b=1595:2783;
                    count_edge_groenland=original(a,b);
                    interG=sum(count_edge_groenland); % somme intermédiaire
                    total_Groenland=sum(interG);   
                   
                    k=1789:2664;
                    l=391:1369;
                    count_edge_Bay=original(k,l);
                    interBay=sum(count_edge_Bay);
                    total_Bay=sum(interBay);   
                   
                    N=1555:2664;
                    n=1600:2762;
                    count_edge_banks=original(N,n);
                    interBanks=sum(count_edge_banks);
                    total_Banks=sum(interBanks);   
                    
                    
                    o=302:1400;
                    p=976:2128;
                    count_edge_labrador=original(o,p);  
                    interLab=sum(count_edge_labrador);
                    total_Labrador=sum(interLab);   
                    
                    % le répertoire de destination DOIT ÊTRE CRÉÉ AVANT l'éxécution 
                    dest ='H:\matlab_atlantic\subset_output\';  
%                         
                     save([dest 'Subset_' year '_' strmonth '.mat'],'total_Groenland','total_Bay','total_Banks','total_Labrador','count_edge_Bay','count_edge_groenland','count_edge_banks','count_edge_labrador')   
%                          save([dest 'Subset_' year '_' strmonth '.mat'],'count_edge_Bay')   

                     clear ('data','count_edge_labrador','count_edge_banks','count_edge_Bay','count_edge_groenland','a','b','k','l','N','n','o','p');
%                         clear ('data','count_edge_Bay','k','l');

                end      
            end    
    end
end       
  
            