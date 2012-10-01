#include <cstring>
#include <seqan.h>
#include <seqan/align.h>
#include <string>

#ifdef SWIG
namespace std
{	
%typemap(in,numinputs=0) std::string& str1(std::string temp) 
	{
		$1 = &temp;
	}
%typemap(in,numinputs=0) std::string& str2(std::string temp) 
	{
		$1 = &temp;
	}
	%typemap(argout)std::string& str1
	{
		$result=SWIG_Python_AppendOutput($result, PyString_FromStringAndSize((*$1).c_str(),(*$1).length()));
	}
	%typemap(argout)std::string& str2
	{
		$result=SWIG_Python_AppendOutput($result, PyString_FromStringAndSize((*$1).c_str(),(*$1).length()));
	}
%typemap(out)std::string
	{
		$result=SWIG_Python_AppendOutput($result, PyString_FromStringAndSize(($1).c_str(),($1).length()));
	}
	
}	


#endif

namespace seqan
{

	#ifdef SWIG

	//typemap(in) -> python to c++
	//typemap(out) -> c++ to python

	//typemap(in,numinputs=0)int& scoreAlignment. This typemap masked scoreAlignment as a shadows argument, which means SWIG will not consider this variable to be filled when calling the function in python. Instead it is filling the argument with a temporary value. The currect value will be replaced in the typemap(argout) later. 
	%typemap(in,numinputs=0)int& scoreAlignment(int temp)
	{
		$1=&temp;
	}
	//this template shows how to return a tupel into python
	%typemap(argout)int& scoreAlignment
	{    
		PyObject *o, *o2, *o3; 
    o = PyInt_FromLong(*$1); //$1 -> the value which has to be convert to python
    if ((!$result) || ($result == Py_None))  //check if already exist a result, if not set the actual result
		{
			$result = o;
    }
		else
		{
      if(!PyTuple_Check($result)) //build a tuple in python 
			 {
            PyObject *o2 = $result;
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
				o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
	}
	
	/*
	//see above... by using char** we have to pre allocate memory
	%typemap(in,numinputs=0) char** str1(char* temp) 
	{
		$1 = (char **) malloc((strlen(temp)+1)*sizeof(char *));
	}
	%typemap(in,numinputs=0) char** str2(char* temp) 
	{
		  $1 = (char **) malloc((strlen(temp)+1)*sizeof(char *));
	}
	//returns a list of char using the Python_AppendOupt function.
	*/

	
	/*
	%typemap(argout)char** str1
	{
		$result=SWIG_Python_AppendOutput($result, PyString_FromStringAndSize(*$1,strlen(*$1)));
	}
	%typemap(argout)char** str2
	{
		$result=SWIG_Python_AppendOutput($result,PyString_FromString(*$1));
	}
	//destroy the pre allocate memory -> avoid memory leaks
	%typemap(freearg) char ** 
	{
		free((char *) $1);
	}
	*/
	%typemap(in,numinputs=0) int** scoreMatrix(int* temp)
	{
		  $1 = (int **) malloc((sizeof(temp)+1)*sizeof(int *));
	}
	%typemap(in,numinputs=0) int& aSize(int temp)
	{
		  $1 = &temp;
	}
	//int array -> list of integer in python. Two information are important, the array itself and its size.
	%typemap(argout) (int** scoreMatrix, int& aSize)
	{
  	int i;
  	$result = PyList_New(*$2); //build a list with the the size aSize
  	for (i = 0; i < *$2; i++)
	 {
    	PyObject *o =PyInt_FromLong(*(*$1+i));
    	PyList_SetItem($result,i,o);
  	}
	}


	%typemap(freearg) int ** 
	{
		free((int *) $1);
	}

	//python integer list -> int array
	%typemap(in) (int len, int *value)
	{
		int i;
		if (!PyList_Check($input))
		{
		  PyErr_SetString(PyExc_ValueError, "Expecting a list");
		  return NULL;
		}
		$1 = PyList_Size($input); //get size of the list
		$2 = (int*) malloc(($1+1)*sizeof(int)); //allocate memory with the correct size
		for (i = 0; i < $1; i++)
		{
		  PyObject *s = PyList_GetItem($input,i);
		  if (!PyInt_Check(s))
			{
		      free($2);
		      PyErr_SetString(PyExc_ValueError, "List items must be integers");
		      return NULL;
		  }
		  $2[i] = (int)PyInt_AS_LONG(s); //set the value into the array
		}
		$2[i] = 0;
	}

	%typemap(freearg) (int len, int *value) 
	{
		 if ($2) free($2);
	}
	//shows how to set default values in SWIG
	%typemap(default) unsigned int matrix
	{
		 $1 = 0;
	}
	%typemap(default) unsigned int dia1
	{
		 $1 = 0;
	}
	%typemap(default) unsigned int dia2 
	{
		 $1 = 0;
	}


	%define DOCSTRING_
	"getAminoAcidScoreMatrix(\"X\")->[long integer]; \n Returns a list containing integer, which represent the specified scoring matrices X.\n Selectable value for X are:\n Blosum30 \n Blosum45 \n Blosum62 \n Blosum80 \n Pam40 \n Pam120 \n Pam200 \n Pam250 \n Vtml200 \n Example: scoreMatrix=getAminoAcidScoreMatrix(\"Blosum30\");\n print(scoreMatrix)"
	%enddef

	%define DOCSTRING2
	"printAlignment(AlignmentObject)->[Strings];\n Returns a list containing three string objects.\n The first string illustrate the alignment.\n String two and three contains the alignment sequence separately.\n"
	%enddef

	%feature("autodoc", DOCSTRING_) getAminoAcidScoreMatrix;

	%feature("autodoc", DOCSTRING2) printAlignment;
	
	#endif
	
	unsigned int binToInt(std::string s)
	{
	  unsigned int lastnumber = s[s.length()-1] == '1' ? 1 : 0;
	  if(s.length() > 1)
	  {
		  return 2 * binToInt(s.substr(0,s.length() - 1)) + lastnumber;
	  }
	  else
	  {
		  return lastnumber;
	  }

	}

	//returns a specific amino acid scoring matrix selected by the input name.
	void getAminoAcidScoreMatrix(char* str,int** scoreMatrix, int& aSize)
	{
		if(std::strcmp(str,"Blosum30")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Blosum30> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Blosum30 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix); 
			aSize=TScore::TAB_SIZE;
		}

		else if (std::strcmp(str,"Blosum45")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Blosum45> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Blosum45 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Blosum62")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Blosum62> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Blosum62 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Blosum80")==0)
		{

			typedef Score<int, ScoreMatrix<AminoAcid, _Blosum80> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Blosum80 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;

		}
		else if(std::strcmp(str,"Pam40")==0)
		{
	
			typedef Score<int, ScoreMatrix<AminoAcid, _Pam40> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Pam40 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Pam120")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Pam120> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Pam120 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Pam200")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Pam200> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Pam200 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Pam250")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Pam250> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Pam250 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
		else if(std::strcmp(str,"Vtml200")==0)
		{
			typedef Score<int, ScoreMatrix<AminoAcid, _Vtml200> > TScore;
			int const* matrix=_ScoringMatrixData<int ,AminoAcid, _Vtml200 >::getData();
			*scoreMatrix = new int[TScore::TAB_SIZE];
			arrayCopy(matrix, matrix + TScore::TAB_SIZE, *scoreMatrix);
			aSize=TScore::TAB_SIZE;
		}
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////

//returns the score of an alignment and save the alignment itself using a reference on AlignObject.
//char algo determines what kind of alignment algorithm is chosen.
	template<typename Type,typename TSpec, typename TAlignConfig>
	int _alignment(seqan::Align<seqan::String<Type> >*align,seqan::Score<int,TSpec>& score,TAlignConfig& ac,unsigned int dia1, unsigned int dia2,char *algo)
	{
		if(std::strcmp(algo,"NeedlemanWunsch")==0) return globalAlignment(*align,score,ac,NeedlemanWunsch());
		else if(std::strcmp(algo,"Gotoh")==0) return globalAlignment(*align,score,ac,Gotoh());
		else if(std::strcmp(algo,"BandedNeedlemanWunsch")==0) return  globalAlignment(stringSet(*align),score,ac,dia1,dia2,BandedNeedlemanWunsch());
		else if(std::strcmp(algo,"BandedGotoh")==0) return globalAlignment(stringSet(*align),score,ac,dia1,dia2,BandedGotoh());
	}
	 
	//returns the score of an alignment by calling the above _alignment() function. This function determines the
	//correct AlignConfig using the integer given by the user.
	template<typename Type, typename TSpec>
	int _alignment(seqan::Align<seqan::String<Type> >*align,seqan::Score<int,TSpec>& score,unsigned int ac,char * algo,unsigned int dia1, unsigned int dia2)
	{
	int alignmentScore=0;
			switch(ac)
			{
				case 0: {
									AlignConfig<false,false,false,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);
									break;
								}
				case 1: {
									AlignConfig<false,false,false,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);
									break;
								}
				case 2: {
									AlignConfig<false,false,true,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);
									break;
								}
				case 3: {
									AlignConfig<false,false,true,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 4: {
									AlignConfig<false,true,false,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 5: {
									AlignConfig<false,true,false,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 6: {
									AlignConfig<false,true,true,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 7: {
									AlignConfig<false,true,true,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 8: {
									AlignConfig<true,false,false,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 9: {
									AlignConfig<true,false,false,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 10: {
									AlignConfig<true,false,true,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 11: {
									AlignConfig<true,false,true,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 12: {
									AlignConfig<true,true,false,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 13: {
									AlignConfig<true,true,false,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 14: {
									AlignConfig<true,true,true,false>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
				case 15: {
									AlignConfig<true,true,true,true>alignconfig;
									alignmentScore=_alignment(align,score,alignconfig,dia1,dia2,algo);							
									break;
								}
			}
	return alignmentScore;
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//main wrapper function for calling the alignment.
	template< typename Type >
	seqan::Align<seqan::String<Type> >* alignSequence(char* seq1, char*seq2, int len, int *value, char* algo, char* binary, unsigned int dia1, unsigned int dia2, int& scoreAlignment)
	{
		unsigned int matrix=0; 
		std::string binarySt=binary;
		if(binarySt.length()>1)
		matrix= binToInt(binarySt.substr(2,binarySt.length()));
		else matrix=0;
		
		if(matrix >15) matrix=0;
		Align<String<Type> >* align = new Align<String<Type> >();
		
		String<Type> _seq1=seq1; //this is necessary!!! Python strings are immutable!! Swig does not care what happens with the strings
		String<Type> _seq2=seq2;
		
		resize(rows(*align), 2);
		assignSource(row(*align, 0), _seq1);
		assignSource(row(*align, 1), _seq2);
		//check if the actual chosen alignment algorithm supports a score matrix
		if(len >=4 && ValueSize<Type>::VALUE * ValueSize<Type>::VALUE == len && std::strcmp(algo,"MyersHirschberg")!=0 && std::strcmp(algo,"Hirschberg")!=0 && std::strcmp(algo,"SmithWaterman")!=0)
		{
			Score<int, ScoreMatrix<Type> > score;
 			arrayCopy(value, value + len, score.data_tab);
			scoreAlignment=_alignment(align,score,matrix,algo,dia1,dia2);
		}
		else if(len == 4) 
		{   
			Score<int, Simple > score(*(value),*(value+1),*(value+2),*(value+3));
			//check if we have to set AlignConfig for the current alignment algorithm
			if(std::strcmp(algo,"MyersHirschberg")!=0 && std::strcmp(algo,"Hirschberg")!=0&&std::strcmp(algo,"SmithWaterman")!=0 )
				{ 
					scoreAlignment=_alignment(align,score,matrix,algo,dia1,dia2);
				}
				else if(std::strcmp(algo,"MyersHirschberg")==0)
				{
					scoreAlignment=globalAlignment(*align,score,MyersHirschberg());
				}
				else if(std::strcmp(algo,"Hirschberg")==0)
				{
					scoreAlignment=globalAlignment(*align,score,Hirschberg());		
				}
				else if(std::strcmp(algo,"SmithWaterman")==0)
				{
					scoreAlignment=localAlignment(*align,score, SmithWaterman());
				}
		}	
//std::cout << *align << std::endl;
		return align;
	} 
////////////////////////////////////////////////////////////////////////////////////////////////////////

	//returns the alignment into different output variables ->via typemap into a string list
	template<typename Type>
	std::string printAlignment(seqan::Align<seqan::String<Type> >*align, std::string& str1, std::string& str2)//char **str1,char **str2)
	{
		typedef Align<seqan::String<Type> > const TAlign;
		typedef typename Row<TAlign>::Type TRow;
		typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
		typedef typename Position<TAlign>::Type TPosition;
		TRowsPosition row_count = length(rows(*align));
		TPosition begin_ = beginPosition(cols(*align));
		TPosition end_ = endPosition(cols(*align));
		unsigned int baseCount=0;
		unsigned int leftSpace=6;
		std::string pAlignment;
		std::string seq1Alignment,seq2Alignment;
		while(begin_ < end_) 
		{
			unsigned int windowSize_ = 50;
			if ((begin_ + windowSize_)>end_) windowSize_=end_ - begin_;
			for(TRowsPosition i=0;i<2*row_count-1;++i)
			{
				for(unsigned int j = 0;j<leftSpace+2;++j)
				{
				pAlignment+=" "; 
				}
				if ((i % 2)==0) 
				{
					TRow& row_ = row(*align, i/2);
					typedef typename Iterator<typename Row<TAlign>::Type const, Standard>::Type TIter;
					TIter begin1_ = iter(row_, begin_);
					TIter end1_ = iter(row_, begin_ + windowSize_);
					for (; begin1_ != end1_; ++begin1_) 
					{
						if (isGap(begin1_))
						{
							pAlignment+="-";
							if(i==0) 
							{
								seq1Alignment+="-";
							}
							else if(i==2)seq2Alignment+="-"; 
						}
						else
						{
							char b=convert<char>(*begin1_);
							char* t=(char*)b;
							pAlignment+=b;
							if(i==0) seq1Alignment+=b;
							else if(i==2)seq2Alignment +=b;
						}
					}
				}
				else
				{
					for(unsigned int j = 0;j<windowSize_;++j)
					{
						if ((!isGap(row(*align, (i-1)/2), begin_+j)) &&(!isGap(row(*align, (i+1)/2), begin_+j)) && (row(*align, (i-1)/2)[begin_+j]==row(*align, (i+1)/2)[begin_+j]))
						{
							pAlignment+="|";
						}
						else
					 {
						pAlignment+=" ";
						}
					} 
				}
				pAlignment+="\n";
			}
			pAlignment+="\n";
			begin_+=50;
		}
	pAlignment+="\n";

	str1= seq1Alignment;
	str2= seq2Alignment;
	return pAlignment;
	}
/*
//overload function with seqan pointers does not work......
	char* printAlignment(seqan::Align<seqan::String<seqan::Dna5> >&align,char **str1,char **str2)
	{
		return _printAlignment(align,str1,str2);
	}

	char* printAlignment(seqan::Align<seqan::String<seqan::AminoAcid> >*align,char **str1,char **str2)
	{
		
		return _printAlignment(*align,str1,str2);
	}

	char* printAlignment(seqan::Align<seqan::String<char> >&align,char **str1,char **str2)
	{
		return _printAlignment(align,str1,str2);
	}
*/
}
