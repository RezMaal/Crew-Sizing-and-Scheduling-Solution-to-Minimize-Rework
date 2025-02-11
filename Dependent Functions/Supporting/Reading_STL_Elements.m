%% Prepared, Written and Created by Jun.-Prof. Dr. Reza Maalek
% Only for registered members of the "Digital Technologies for Field
% Information Modeling" Course at KIT. Use of the software for commercial
% purposes is strictly prohibited. Distribution only allowed with written
% concent from Jun.-Prof. Dr. Reza Maalek (reza.maalek@kit.edu). For
% educational use only.

function stz=Reading_STL_Elements(maxNum,OBJ)
    str="ObjX.stl";
    L=length(OBJ);
    stz=zeros(sum(maxNum),1);
    Mn=[0,maxNum(1:end-1)];
    stz=string(stz);
    for i=1:L
        strn=replace(str,'Obj',OBJ(i));
        for Num=1:maxNum(i)
            stz(sum(Mn(1:i))+Num)=replace(strn,'X',string(Num));
        end
    end  
end