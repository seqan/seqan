// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// FASTA-Tests for seqan/stream/read_fasta_fastq.h
// ==========================================================================


#ifndef TEST_STREAM_TEST_STREAM_READ_GENBANK_H_
#define TEST_STREAM_TEST_STREAM_READ_GENBANK_H_

// Return an GenBank format string with two entries.
char const * testHelperReturnGenBankFile()
{
    char const * result =
            "LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999\n"
            "DEFINITION  Saccharomyces cerevisiae TCP1-beta gene, partial cds, and Axl2p\n"
            "            (AXL2) and Rev7p (REV7) genes, complete cds.\n"
            "ACCESSION   U49845\n"
            "VERSION     U49845.1  GI:1293613\n"
            "KEYWORDS    .\n"
            "SOURCE      Saccharomyces cerevisiae (baker's yeast)\n"
            "  ORGANISM  Saccharomyces cerevisiae\n"
            "            Eukaryota; Fungi; Ascomycota; Saccharomycotina; Saccharomycetes;\n"
            "            Saccharomycetales; Saccharomycetaceae; Saccharomyces.\n"
            "REFERENCE   1  (bases 1 to 5028)\n"
            "  AUTHORS   Torpey,L.E., Gibbs,P.E., Nelson,J. and Lawrence,C.W.\n"
            "  TITLE     Cloning and sequence of REV7, a gene whose function is required for\n"
            "            DNA damage-induced mutagenesis in Saccharomyces cerevisiae\n"
            "  JOURNAL   Yeast 10 (11), 1503-1509 (1994)\n"
            "  PUBMED    7871890\n"
            "REFERENCE   2  (bases 1 to 5028)\n"
            "  AUTHORS   Roemer,T., Madden,K., Chang,J. and Snyder,M.\n"
            "  TITLE     Selection of axial growth sites in yeast requires Axl2p, a novel\n"
            "            plasma membrane glycoprotein\n"
            "  JOURNAL   Genes Dev. 10 (7), 777-793 (1996)\n"
            "  PUBMED    8846915\n"
            "REFERENCE   3  (bases 1 to 5028)\n"
            "  AUTHORS   Roemer,T.\n"
            "  TITLE     Direct Submission\n"
            "  JOURNAL   Submitted (22-FEB-1996) Terry Roemer, Biology, Yale University, New\n"
            "            Haven, CT, USA\n"
            "FEATURES             Location/Qualifiers\n"
            "     source          1..5028\n"
            "                     /organism=\"Saccharomyces cerevisiae\"\n"
            "                     /db_xref=\"taxon:4932\"\n"
            "                     /chromosome=\"IX\"\n"
            "                     /map=\"9\"\n"
            "     CDS             <1..206\n"
            "                     /codon_start=3\n"
            "                     /product=\"TCP1-beta\"\n"
            "                     /protein_id=\"AAA98665.1\"\n"
            "                     /db_xref=\"GI:1293614\"\n"
            "                     /translation=\"SSIYNGISTSGLDLNNGTIADMRQLGIVESYKLKRAVVSSASEA\n"
            "                     AEVLLRVDNIIRARPRTANRQHM\"\n"
            "     gene            687..3158\n"
            "                     /gene=\"AXL2\"\n"
            "     CDS             687..3158\n"
            "                     /gene=\"AXL2\"\n"
            "                     /note=\"plasma membrane glycoprotein\"\n"
            "                     /codon_start=1\n"
            "                     /function=\"required for axial budding pattern of S.\n"
            "                     cerevisiae\"\n"
            "                     /product=\"Axl2p\"\n"
            "                     /protein_id=\"AAA98666.1\"\n"
            "                     /db_xref=\"GI:1293615\"\n"
            "                     /translation=\"MTQLQISLLLTATISLLHLVVATPYEAYPIGKQYPPVARVNESF\n"
            "                     TFQISNDTYKSSVDKTAQITYNCFDLPSWLSFDSSSRTFSGEPSSDLLSDANTTLYFN\n"
            "                     VILEGTDSADSTSLNNTYQFVVTNRPSISLSSDFNLLALLKNYGYTNGKNALKLDPNE\n"
            "                     VFNVTFDRSMFTNEESIVSYYGRSQLYNAPLPNWLFFDSGELKFTGTAPVINSAIAPE\n"
            "                     TSYSFVIIATDIEGFSAVEVEFELVIGAHQLTTSIQNSLIINVTDTGNVSYDLPLNYV\n"
            "                     YLDDDPISSDKLGSINLLDAPDWVALDNATISGSVPDELLGKNSNPANFSVSIYDTYG\n"
            "                     DVIYFNFEVVSTTDLFAISSLPNINATRGEWFSYYFLPSQFTDYVNTNVSLEFTNSSQ\n"
            "                     DHDWVKFQSSNLTLAGEVPKNFDKLSLGLKANQGSQSQELYFNIIGMDSKITHSNHSA\n"
            "                     NATSTRSSHHSTSTSSYTSSTYTAKISSTSAAATSSAPAALPAANKTSSHNKKAVAIA\n"
            "                     CGVAIPLGVILVALICFLIFWRRRRENPDDENLPHAISGPDLNNPANKPNQENATPLN\n"
            "                     NPFDDDASSYDDTSIARRLAALNTLKLDNHSATESDISSVDEKRDSLSGMNTYNDQFQ\n"
            "                     SQSKEELLAKPPVQPPESPFFDPQNRSSSVYMDSEPAVNKSWRYTGNLSPVSDIVRDS\n"
            "                     YGSQKTVDTEKLFDLEAPEKEKRTSRDVTMSSLDPWNSNISPSPVRKSVTPSPYNVTK\n"
            "                     HRNRHLQNIQDSQSGKNGITPTTMSTSSSDDFVPVKDGENFCWVHSMEPDRRPSKKRL\n"
            "                     VDFSNKSNVNVGQVKDIHGRIPEML\"\n"
            "     gene            complement(3300..4037)\n"
            "                     /gene=\"REV7\"\n"
            "     CDS             complement(3300..4037)\n"
            "                     /gene=\"REV7\"\n"
            "                     /codon_start=1\n"
            "                     /product=\"Rev7p\"\n"
            "                     /protein_id=\"AAA98667.1\"\n"
            "                     /db_xref=\"GI:1293616\"\n"
            "                     /translation=\"MNRWVEKWLRVYLKCYINLILFYRNVYPPQSFDYTTYQSFNLPQ\n"
            "                     FVPINRHPALIDYIEELILDVLSKLTHVYRFSICIINKKNDLCIEKYVLDFSELQHVD\n"
            "                     KDDQIITETEVFDEFRSSLNSLIMHLEKLPKVNDDTITFEAVINAIELELGHKLDRNR\n"
            "                     RVDSLEEKAEIERDSNWVKCQEDENLPDNNGFQPPKIKLTSLVGSDVGPLIIHQFSEK\n"
            "                     LISGDDKILNGVYSQYEEGESIFGSLF\"\n"
            "ORIGIN\n"
            "        1 gatcctccat atacaacggt atctccacct caggtttaga tctcaacaac ggaaccattg\n"
            "       61 ccgacatgag acagttaggt atcgtcgaga gttacaagct aaaacgagca gtagtcagct\n"
            "      121 ctgcatctga agccgctgaa gttctactaa gggtggataa catcatccgt gcaagaccaa\n"
            "      181 gaaccgccaa tagacaacat atgtaacata tttaggatat acctcgaaaa taataaaccg\n"
            "      241 ccacactgtc attattataa ttagaaacag aacgcaaaaa ttatccacta tataattcaa\n"
            "      301 agacgcgaaa aaaaaagaac aacgcgtcat agaacttttg gcaattcgcg tcacaaataa\n"
            "      361 attttggcaa cttatgtttc ctcttcgagc agtactcgag ccctgtctca agaatgtaat\n"
            "      421 aatacccatc gtaggtatgg ttaaagatag catctccaca acctcaaagc tccttgccga\n"
            "      481 gagtcgccct cctttgtcga gtaattttca cttttcatat gagaacttat tttcttattc\n"
            "      541 tttactctca catcctgtag tgattgacac tgcaacagcc accatcacta gaagaacaga\n"
            "      601 acaattactt aatagaaaaa ttatatcttc ctcgaaacga tttcctgctt ccaacatcta\n"
            "      661 cgtatatcaa gaagcattca cttaccatga cacagcttca gatttcatta ttgctgacag\n"
            "      721 ctactatatc actactccat ctagtagtgg ccacgcccta tgaggcatat cctatcggaa\n"
            "      781 aacaataccc cccagtggca agagtcaatg aatcgtttac atttcaaatt tccaatgata\n"
            "      841 cctataaatc gtctgtagac aagacagctc aaataacata caattgcttc gacttaccga\n"
            "      901 gctggctttc gtttgactct agttctagaa cgttctcagg tgaaccttct tctgacttac\n"
            "      961 tatctgatgc gaacaccacg ttgtatttca atgtaatact cgagggtacg gactctgccg\n"
            "     1021 acagcacgtc tttgaacaat acataccaat ttgttgttac aaaccgtcca tccatctcgc\n"
            "     1081 tatcgtcaga tttcaatcta ttggcgttgt taaaaaacta tggttatact aacggcaaaa\n"
            "     1141 acgctctgaa actagatcct aatgaagtct tcaacgtgac ttttgaccgt tcaatgttca\n"
            "     1201 ctaacgaaga atccattgtg tcgtattacg gacgttctca gttgtataat gcgccgttac\n"
            "     1261 ccaattggct gttcttcgat tctggcgagt tgaagtttac tgggacggca ccggtgataa\n"
            "     1321 actcggcgat tgctccagaa acaagctaca gttttgtcat catcgctaca gacattgaag\n"
            "     1381 gattttctgc cgttgaggta gaattcgaat tagtcatcgg ggctcaccag ttaactacct\n"
            "     1441 ctattcaaaa tagtttgata atcaacgtta ctgacacagg taacgtttca tatgacttac\n"
            "     1501 ctctaaacta tgtttatctc gatgacgatc ctatttcttc tgataaattg ggttctataa\n"
            "     1561 acttattgga tgctccagac tgggtggcat tagataatgc taccatttcc gggtctgtcc\n"
            "     1621 cagatgaatt actcggtaag aactccaatc ctgccaattt ttctgtgtcc atttatgata\n"
            "     1681 cttatggtga tgtgatttat ttcaacttcg aagttgtctc cacaacggat ttgtttgcca\n"
            "     1741 ttagttctct tcccaatatt aacgctacaa ggggtgaatg gttctcctac tattttttgc\n"
            "     1801 cttctcagtt tacagactac gtgaatacaa acgtttcatt agagtttact aattcaagcc\n"
            "     1861 aagaccatga ctgggtgaaa ttccaatcat ctaatttaac attagctgga gaagtgccca\n"
            "     1921 agaatttcga caagctttca ttaggtttga aagcgaacca aggttcacaa tctcaagagc\n"
            "     1981 tatattttaa catcattggc atggattcaa agataactca ctcaaaccac agtgcgaatg\n"
            "     2041 caacgtccac aagaagttct caccactcca cctcaacaag ttcttacaca tcttctactt\n"
            "     2101 acactgcaaa aatttcttct acctccgctg ctgctacttc ttctgctcca gcagcgctgc\n"
            "     2161 cagcagccaa taaaacttca tctcacaata aaaaagcagt agcaattgcg tgcggtgttg\n"
            "     2221 ctatcccatt aggcgttatc ctagtagctc tcatttgctt cctaatattc tggagacgca\n"
            "     2281 gaagggaaaa tccagacgat gaaaacttac cgcatgctat tagtggacct gatttgaata\n"
            "     2341 atcctgcaaa taaaccaaat caagaaaacg ctacaccttt gaacaacccc tttgatgatg\n"
            "     2401 atgcttcctc gtacgatgat acttcaatag caagaagatt ggctgctttg aacactttga\n"
            "     2461 aattggataa ccactctgcc actgaatctg atatttccag cgtggatgaa aagagagatt\n"
            "     2521 ctctatcagg tatgaataca tacaatgatc agttccaatc ccaaagtaaa gaagaattat\n"
            "     2581 tagcaaaacc cccagtacag cctccagaga gcccgttctt tgacccacag aataggtctt\n"
            "     2641 cttctgtgta tatggatagt gaaccagcag taaataaatc ctggcgatat actggcaacc\n"
            "     2701 tgtcaccagt ctctgatatt gtcagagaca gttacggatc acaaaaaact gttgatacag\n"
            "     2761 aaaaactttt cgatttagaa gcaccagaga aggaaaaacg tacgtcaagg gatgtcacta\n"
            "     2821 tgtcttcact ggacccttgg aacagcaata ttagcccttc tcccgtaaga aaatcagtaa\n"
            "     2881 caccatcacc atataacgta acgaagcatc gtaaccgcca cttacaaaat attcaagact\n"
            "     2941 ctcaaagcgg taaaaacgga atcactccca caacaatgtc aacttcatct tctgacgatt\n"
            "     3001 ttgttccggt taaagatggt gaaaattttt gctgggtcca tagcatggaa ccagacagaa\n"
            "     3061 gaccaagtaa gaaaaggtta gtagattttt caaataagag taatgtcaat gttggtcaag\n"
            "     3121 ttaaggacat tcacggacgc atcccagaaa tgctgtgatt atacgcaacg atattttgct\n"
            "     3181 taattttatt ttcctgtttt attttttatt agtggtttac agatacccta tattttattt\n"
            "     3241 agtttttata cttagagaca tttaatttta attccattct tcaaatttca tttttgcact\n"
            "     3301 taaaacaaag atccaaaaat gctctcgccc tcttcatatt gagaatacac tccattcaaa\n"
            "     3361 attttgtcgt caccgctgat taatttttca ctaaactgat gaataatcaa aggccccacg\n"
            "     3421 tcagaaccga ctaaagaagt gagttttatt ttaggaggtt gaaaaccatt attgtctggt\n"
            "     3481 aaattttcat cttcttgaca tttaacccag tttgaatccc tttcaatttc tgctttttcc\n"
            "     3541 tccaaactat cgaccctcct gtttctgtcc aacttatgtc ctagttccaa ttcgatcgca\n"
            "     3601 ttaataactg cttcaaatgt tattgtgtca tcgttgactt taggtaattt ctccaaatgc\n"
            "     3661 ataatcaaac tatttaagga agatcggaat tcgtcgaaca cttcagtttc cgtaatgatc\n"
            "     3721 tgatcgtctt tatccacatg ttgtaattca ctaaaatcta aaacgtattt ttcaatgcat\n"
            "     3781 aaatcgttct ttttattaat aatgcagatg gaaaatctgt aaacgtgcgt taatttagaa\n"
            "     3841 agaacatcca gtataagttc ttctatatag tcaattaaag caggatgcct attaatggga\n"
            "     3901 acgaactgcg gcaagttgaa tgactggtaa gtagtgtagt cgaatgactg aggtgggtat\n"
            "     3961 acatttctat aaaataaaat caaattaatg tagcatttta agtataccct cagccacttc\n"
            "     4021 tctacccatc tattcataaa gctgacgcaa cgattactat tttttttttc ttcttggatc\n"
            "     4081 tcagtcgtcg caaaaacgta taccttcttt ttccgacctt ttttttagct ttctggaaaa\n"
            "     4141 gtttatatta gttaaacagg gtctagtctt agtgtgaaag ctagtggttt cgattgactg\n"
            "     4201 atattaagaa agtggaaatt aaattagtag tgtagacgta tatgcatatg tatttctcgc\n"
            "     4261 ctgtttatgt ttctacgtac ttttgattta tagcaagggg aaaagaaata catactattt\n"
            "     4321 tttggtaaag gtgaaagcat aatgtaaaag ctagaataaa atggacgaaa taaagagagg\n"
            "     4381 cttagttcat cttttttcca aaaagcaccc aatgataata actaaaatga aaaggatttg\n"
            "     4441 ccatctgtca gcaacatcag ttgtgtgagc aataataaaa tcatcacctc cgttgccttt\n"
            "     4501 agcgcgtttg tcgtttgtat cttccgtaat tttagtctta tcaatgggaa tcataaattt\n"
            "     4561 tccaatgaat tagcaatttc gtccaattct ttttgagctt cttcatattt gctttggaat\n"
            "     4621 tcttcgcact tcttttccca ttcatctctt tcttcttcca aagcaacgat ccttctaccc\n"
            "     4681 atttgctcag agttcaaatc ggcctctttc agtttatcca ttgcttcctt cagtttggct\n"
            "     4741 tcactgtctt ctagctgttg ttctagatcc tggtttttct tggtgtagtt ctcattatta\n"
            "     4801 gatctcaagt tattggagtc ttcagccaat tgctttgtat cagacaattg actctctaac\n"
            "     4861 ttctccactt cactgtcgag ttgctcgttt ttagcggaca aagatttaat ctcgttttct\n"
            "     4921 ttttcagtgt tagattgctc taattctttg agctgttctc tcagctcctc atatttttct\n"
            "     4981 tgccatgact cagattctaa ttttaagcta ttcaatttct ctttgatc\n"
            "//\n"
            "LOCUS       AB031069     2487 bp    mRNA            PRI       27-MAY-2000\n"
            "DEFINITION  Homo sapiens PCCX1 mRNA for protein containing CXXC domain 1,\n"
            "            complete cds.\n"
            "ACCESSION   AB031069\n"
            "VERSION     AB031069.1  GI:8100074\n"
            "KEYWORDS    .\n"
            "SOURCE      Homo sapiens embryo male lung fibroblast cell_line:HuS-L12 cDNA to\n"
            "            mRNA.\n"
            "  ORGANISM  Homo sapiens\n"
            "            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;\n"
            "            Mammalia; Eutheria; Primates; Catarrhini; Hominidae; Homo.\n"
            "REFERENCE   1  (sites)\n"
            "  AUTHORS   Fujino,T., Hasegawa,M., Shibata,S., Kishimoto,T., Imai,Si. and\n"
            "            Takano,T.\n"
            "  TITLE     PCCX1, a novel DNA-binding protein with PHD finger and CXXC domain,\n"
            "            is regulated by proteolysis\n"
            "  JOURNAL   Biochem. Biophys. Res. Commun. 271 (2), 305-310 (2000)\n"
            "  MEDLINE   20261256\n"
            "REFERENCE   2  (bases 1 to 2487)\n"
            "  AUTHORS   Fujino,T., Hasegawa,M., Shibata,S., Kishimoto,T., Imai,S. and\n"
            "            Takano,T.\n"
            "  TITLE     Direct Submission\n"
            "  JOURNAL   Submitted (15-AUG-1999) to the DDBJ/EMBL/GenBank databases.\n"
            "            Tadahiro Fujino, Keio University School of Medicine, Department of\n"
            "            Microbiology; Shinanomachi 35, Shinjuku-ku, Tokyo 160-8582, Japan\n"
            "            (E-mail:fujino@microb.med.keio.ac.jp,\n"
            "            Tel:+81-3-3353-1211(ex.62692), Fax:+81-3-5360-1508)\n"
            "FEATURES             Location/Qualifiers\n"
            "     source          1..2487\n"
            "                     /organism=\"Homo sapiens\"\n"
            "                     /db_xref=\"taxon:9606\"\n"
            "                     /sex=\"male\"\n"
            "                     /cell_line=\"HuS-L12\"\n"
            "                     /cell_type=\"lung fibroblast\"\n"
            "                     /dev_stage=\"embryo\"\n"
            "     gene            229..2199\n"
            "                     /gene=\"PCCX1\"\n"
            "     CDS             229..2199\n"
            "                     /gene=\"PCCX1\"\n"
            "                     /note=\"a nuclear protein carrying a PHD finger and a CXXC\n"
            "                     domain\"\n"
            "                     /codon_start=1\n"
            "                     /product=\"protein containing CXXC domain 1\"\n"
            "                     /protein_id=\"BAA96307.1\"\n"
            "                     /db_xref=\"GI:8100075\"\n"
            "                     /translation=\"MEGDGSDPEPPDAGEDSKSENGENAPIYCICRKPDINCFMIGCD\n"
            "                     NCNEWFHGDCIRITEKMAKAIREWYCRECREKDPKLEIRYRHKKSRERDGNERDSSEP\n"
            "                     RDEGGGRKRPVPDPDLQRRAGSGTGVGAMLARGSASPHKSSPQPLVATPSQHHQQQQQ\n"
            "                     QIKRSARMCGECEACRRTEDCGHCDFCRDMKKFGGPNKIRQKCRLRQCQLRARESYKY\n"
            "                     FPSSLSPVTPSESLPRPRRPLPTQQQPQPSQKLGRIREDEGAVASSTVKEPPEATATP\n"
            "                     EPLSDEDLPLDPDLYQDFCAGAFDDHGLPWMSDTEESPFLDPALRKRAVKVKHVKRRE\n"
            "                     KKSEKKKEERYKRHRQKQKHKDKWKHPERADAKDPASLPQCLGPGCVRPAQPSSKYCS\n"
            "                     DDCGMKLAANRIYEILPQRIQQWQQSPCIAEEHGKKLLERIRREQQSARTRLQEMERR\n"
            "                     FHELEAIILRAKQQAVREDEESNEGDSDDTDLQIFCVSCGHPINPRVALRHMERCYAK\n"
            "                     YESQTSFGSMYPTRIEGATRLFCDVYNPQSKTYCKRLQVLCPEHSRDPKVPADEVCGC\n"
            "                     PLVRDVFELTGDFCRLPKRQCNRHYCWEKLRRAEVDLERVRVWYKLDELFEQERNVRT\n"
            "                     AMTNRAGLLALMLHQTIQHDPLTTDLRSSADR\"\n"
            "BASE COUNT      564 a    715 c    768 g    440 t\n"
            "ORIGIN      \n"
            "        1 agatggcggc gctgaggggt cttgggggct ctaggccggc cacctactgg tttgcagcgg\n"
            "       61 agacgacgca tggggcctgc gcaataggag tacgctgcct gggaggcgtg actagaagcg\n"
            "      121 gaagtagttg tgggcgcctt tgcaaccgcc tgggacgccg ccgagtggtc tgtgcaggtt\n"
            "      181 cgcgggtcgc tggcgggggt cgtgagggag tgcgccggga gcggagatat ggagggagat\n"
            "      241 ggttcagacc cagagcctcc agatgccggg gaggacagca agtccgagaa tggggagaat\n"
            "      301 gcgcccatct actgcatctg ccgcaaaccg gacatcaact gcttcatgat cgggtgtgac\n"
            "      361 aactgcaatg agtggttcca tggggactgc atccggatca ctgagaagat ggccaaggcc\n"
            "      421 atccgggagt ggtactgtcg ggagtgcaga gagaaagacc ccaagctaga gattcgctat\n"
            "      481 cggcacaaga agtcacggga gcgggatggc aatgagcggg acagcagtga gccccgggat\n"
            "      541 gagggtggag ggcgcaagag gcctgtccct gatccagacc tgcagcgccg ggcagggtca\n"
            "      601 gggacagggg ttggggccat gcttgctcgg ggctctgctt cgccccacaa atcctctccg\n"
            "      661 cagcccttgg tggccacacc cagccagcat caccagcagc agcagcagca gatcaaacgg\n"
            "      721 tcagcccgca tgtgtggtga gtgtgaggca tgtcggcgca ctgaggactg tggtcactgt\n"
            "      781 gatttctgtc gggacatgaa gaagttcggg ggccccaaca agatccggca gaagtgccgg\n"
            "      841 ctgcgccagt gccagctgcg ggcccgggaa tcgtacaagt acttcccttc ctcgctctca\n"
            "      901 ccagtgacgc cctcagagtc cctgccaagg ccccgccggc cactgcccac ccaacagcag\n"
            "      961 ccacagccat cacagaagtt agggcgcatc cgtgaagatg agggggcagt ggcgtcatca\n"
            "     1021 acagtcaagg agcctcctga ggctacagcc acacctgagc cactctcaga tgaggaccta\n"
            "     1081 cctctggatc ctgacctgta tcaggacttc tgtgcagggg cctttgatga ccatggcctg\n"
            "     1141 ccctggatga gcgacacaga agagtcccca ttcctggacc ccgcgctgcg gaagagggca\n"
            "     1201 gtgaaagtga agcatgtgaa gcgtcgggag aagaagtctg agaagaagaa ggaggagcga\n"
            "     1261 tacaagcggc atcggcagaa gcagaagcac aaggataaat ggaaacaccc agagagggct\n"
            "     1321 gatgccaagg accctgcgtc actgccccag tgcctggggc ccggctgtgt gcgccccgcc\n"
            "     1381 cagcccagct ccaagtattg ctcagatgac tgtggcatga agctggcagc caaccgcatc\n"
            "     1441 tacgagatcc tcccccagcg catccagcag tggcagcaga gcccttgcat tgctgaagag\n"
            "     1501 cacggcaaga agctgctcga acgcattcgc cgagagcagc agagtgcccg cactcgcctt\n"
            "     1561 caggaaatgg aacgccgatt ccatgagctt gaggccatca ttctacgtgc caagcagcag\n"
            "     1621 gctgtgcgcg aggatgagga gagcaacgag ggtgacagtg atgacacaga cctgcagatc\n"
            "     1681 ttctgtgttt cctgtgggca ccccatcaac ccacgtgttg ccttgcgcca catggagcgc\n"
            "     1741 tgctacgcca agtatgagag ccagacgtcc tttgggtcca tgtaccccac acgcattgaa\n"
            "     1801 ggggccacac gactcttctg tgatgtgtat aatcctcaga gcaaaacata ctgtaagcgg\n"
            "     1861 ctccaggtgc tgtgccccga gcactcacgg gaccccaaag tgccagctga cgaggtatgc\n"
            "     1921 gggtgccccc ttgtacgtga tgtctttgag ctcacgggtg acttctgccg cctgcccaag\n"
            "     1981 cgccagtgca atcgccatta ctgctgggag aagctgcggc gtgcggaagt ggacttggag\n"
            "     2041 cgcgtgcgtg tgtggtacaa gctggacgag ctgtttgagc aggagcgcaa tgtgcgcaca\n"
            "     2101 gccatgacaa accgcgcggg attgctggcc ctgatgctgc accagacgat ccagcacgat\n"
            "     2161 cccctcacta ccgacctgcg ctccagtgcc gaccgctgag cctcctggcc cggacccctt\n"
            "     2221 acaccctgca ttccagatgg gggagccgcc cggtgcccgt gtgtccgttc ctccactcat\n"
            "     2281 ctgtttctcc ggttctccct gtgcccatcc accggttgac cgcccatctg cctttatcag\n"
            "     2341 agggactgtc cccgtcgaca tgttcagtgc ctggtggggc tgcggagtcc actcatcctt\n"
            "     2401 gcctcctctc cctgggtttt gttaataaaa ttttgaagaa accaaaaaaa aaaaaaaaaa\n"
            "     2461 aaaaaaaaaa aaaaaaaaaa aaaaaaa\n"
            "//\n";
    return result;
}

template <typename TRecordReader>
void testHelperReadGenBankSingle(TRecordReader & reader)
{
    using namespace seqan;

    CharString buffer, k, v;

    // First record.

    // Read GENBANK headers, split ID header, skip rest.
    readRecord(k, v, reader, GenBankHeader());
    SEQAN_ASSERT_EQ(k, CharString("LOCUS"));
    SEQAN_ASSERT_EQ(v, CharString("SCU49845     5028 bp    DNA             PLN       21-JUN-1999"));
    while (nextIs(reader, GenBankHeader()))
        readRecord(k, v, reader, GenBankHeader());

    // Read GENBANK sequence.
    SEQAN_ASSERT(nextIs(reader, GenBankSequence()));
    DnaString seq;
    readRecord(seq, reader, GenBankSequence());
    SEQAN_ASSERT_EQ(length(seq), 5028u);

    // Second record.

    // Read GENBANK headers, split ID header, skip rest.
    readRecord(k, v, reader, GenBankHeader());
    SEQAN_ASSERT_EQ(k, CharString("LOCUS"));
    SEQAN_ASSERT_EQ(v, CharString("AB031069     2487 bp    mRNA            PRI       27-MAY-2000"));
    while (nextIs(reader, GenBankHeader()))
        readRecord(k, v, reader, GenBankHeader());

    // Read GENBANK sequence.
    SEQAN_ASSERT(nextIs(reader, GenBankSequence()));
    readRecord(seq, reader, GenBankSequence());
    SEQAN_ASSERT_EQ(length(seq), 2487u);

    SEQAN_ASSERT(atEnd(reader));
}

template <typename TRecordReader>
void testHelperReadGenBankRecord(TRecordReader & reader)
{
    using namespace seqan;

    CharString id;
    Dna5String seq;

    readRecord(id, seq, reader, GenBank());
    SEQAN_ASSERT_EQ(id, CharString("U49845.1  GI:1293613"));
    SEQAN_ASSERT_EQ(length(seq), 5028u);

    readRecord(id, seq, reader, GenBank());
    SEQAN_ASSERT_EQ(id, CharString("AB031069.1  GI:8100074"));
    SEQAN_ASSERT_EQ(length(seq), 2487u);
}

SEQAN_DEFINE_TEST(test_stream_read_genbank_single_char_array_stream)
{
    using namespace seqan;

    typedef CharString TStream;

    char const * genbankString = testHelperReturnGenBankFile();
    TStream stream = genbankString;
    DirectionIterator<TStream, Input>::Type reader = directionIterator(stream, Input());

    testHelperReadGenBankSingle(reader);
}

SEQAN_DEFINE_TEST(test_stream_read_genbank_record_char_array_stream)
{
    using namespace seqan;

    typedef CharString TStream;

    char const * genbankString = testHelperReturnGenBankFile();
    TStream stream = genbankString;
    DirectionIterator<TStream, Input>::Type reader = directionIterator(stream, Input());

    testHelperReadGenBankRecord(reader);
}

SEQAN_DEFINE_TEST(test_stream_read_genbank_single_mmap)
{
    using namespace seqan;

    typedef String<char, MMap<> > TString;

    TString mmapString;
    mmapString = testHelperReturnGenBankFile();
    DirectionIterator<TString, Input>::Type reader = directionIterator(mmapString, Input());

    testHelperReadGenBankSingle(reader);
}

SEQAN_DEFINE_TEST(test_stream_read_genbank_single_batch_mmap)
{
    using namespace seqan;

    typedef String<char, MMap<> > TString;

    TString mmapString;
    mmapString = testHelperReturnGenBankFile();
    DirectionIterator<TString, Input>::Type reader = directionIterator(mmapString, Input());

    testHelperReadGenBankRecord(reader);
}

#endif  // TEST_STREAM_TEST_STREAM_READ_GENBANK_H_
