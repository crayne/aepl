


function condit = MakeMedianCol(condit) 

    condit = AeplUtil.MakeConditMat(condit);
    condit.mat(condit.mat>50)=nan;
    condit.medianCol =  nanmedian(condit.mat,2);

end


