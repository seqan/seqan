import alignment as align


dual=bin(2)
alignObject,alignScore=align.alignDna("ACGTC","AACGTCCC",[3,1,-1,0],"NeedlemanWunsch",dual);
print(alignScore);

dnaList=align.printDnaAlignment(alignObject);
print(dnaList[0]);
print(dnaList[1]);
print(dnaList[2]);


scoreMatrix=align.getAminoAcidScoreMatrix("Blosum30");

alignObject,alignScore=align.alignAmino("LYDVAEYAGVSYQTVSRVV","LYDVEEEGVSYQTQTQTVSRVV",scoreMatrix,"NeedlemanWunsch",dual);

aminoList=align.printAminoAcidAlignment(alignObject);
print(alignScore)
print(aminoList[0]);
print(aminoList[1]);
print(aminoList[2]);


