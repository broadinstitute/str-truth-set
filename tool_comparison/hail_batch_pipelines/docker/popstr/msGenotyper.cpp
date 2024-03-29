#include <iostream>
#include <set>
#include <map>
#include <string>
#include <ctime>
#include <math.h>
#include <vector>
#include <numeric>
#include <sstream>
#include <fstream>
#include <ios>
#include <assert.h>
#include <sys/stat.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/basic.h>
#include <seqan/vcf_io.h>
#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <liblinear.hpp>


namespace msGenotyper
{

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

namespace io = boost::iostreams;
using namespace std;
using namespace seqan;

//structure to store read information
struct AttributeLine {
    string PnId;
    float numOfRepeats;
    float ratioBf;
    float ratioAf;
    unsigned locationShift;
    unsigned mateEditDist;
    float purity;
    float ratioOver20In;
    int label;
    double pValue;
} ;

//For storing marker information
struct Marker {
    string chrom;
    int start;
    int end;
    string motif;
    float refRepeatNum;
} ;

//For storing number of members in each class and their pValue-sum
struct LabelProps {
    Pair<int, double> p1;
    Pair<int, double> p2;
    Pair<int, double> p3;
} ;

//For storing marker values
struct MarkerStats
{
    model* regressionModel;
    LabelProps slippCount;
    double slippage;
    unsigned nAlleles;
    double stutter;
    double stepSum;
    std::set<float> alleles;
    std::set<float> trueAlleles;
    double posSlippSum;
    double negSlippSum;
};

//For storing all possible genotypes, their pValues, the chosen genotype, its pValue, all alleles with their frequencies and number of reads available for the decision
struct GenotypeInfo {
    String<Pair<float> > genotypes;
    String<long double> pValues;
    Pair<float> genotype;
    double pValue;
    int numOfReads;
    map<float, int> alleleToFreq; //This maps from reported alleles to their frequencies
    double pValueSum;
    double fullMotifSlippageSum;
} ;

struct MakeGenotypesRet {
	std::set<Pair<float> > genotypesSet;
	String<Pair<float> > genotypes;
};

//So I can map from Markers
bool operator<(const Marker & left, const Marker & right)
{
    return left.start < right.start;
}

//stores various marker specific values
map<Marker, MarkerStats> markerToStats;
//for debugging purposes
bool useGeom = true;

//Parameter, problem and model structs to use in training of logistic regression model and computing pValues
parameter param;
problem prob;
problem probBig; //problem structure including reads with label = 2 to use for predict function
model* model_;
double bias = -1;

//Fills in the x-part of a problem structure from an AttributeLine structure
void fillProblemX(int idx, AttributeLine currentLine, problem& myProb)
{
    myProb.x[idx][0].index = 1;
    myProb.x[idx][0].value = currentLine.ratioBf;
    myProb.x[idx][1].index = 2;
    myProb.x[idx][1].value = currentLine.ratioAf;
    myProb.x[idx][2].index = 3;
    myProb.x[idx][2].value = currentLine.locationShift;
    myProb.x[idx][3].index = 4;
    myProb.x[idx][3].value = currentLine.mateEditDist;
    myProb.x[idx][4].index = 5;
    myProb.x[idx][4].value = currentLine.purity;
    myProb.x[idx][5].index = 6;
    myProb.x[idx][5].value = currentLine.ratioOver20In;
    myProb.x[idx][6].index = -1; // This is to indicate that there aren't any more attributes to read in.
    myProb.x[idx][6].value = 0;
}

double getPval(Marker marker, AttributeLine currentLine)
{
    double predict_label;
    model* model_ = markerToStats[marker].regressionModel;
    double *prob_estimates = (double *) malloc(2*sizeof(double));
    problem prob_;
    prob_.bias = -1;
    prob_.l = 1;
    prob_.n = 6;
    prob_.x = (feature_node **) malloc(prob_.l * sizeof(feature_node *));
    prob_.x[0] = (feature_node *) malloc(10 * sizeof(feature_node));
    fillProblemX(0, currentLine, prob_);
    predict_label = predict_probability(model_,prob_.x[0],prob_estimates);
    return prob_estimates[0];
}

//Parses one line from attribute file by filling up and returning an AttributeLine, also initializes markerToStats map labelProps using the labels
AttributeLine parseNextLine(float winner, float second, bool is_gz, ifstream& attributeFile, io::filtering_istream& attributeFile_gz, Marker& marker, string PnId, map<Pair<string,Marker>, GenotypeInfo>& PnAndMarkerToGenotype, String<string> firstLine, bool useFirstLine, bool useModelAndLabels, bool enoughReads)
{
    PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].genotype = Pair<float>(winner,second);
    AttributeLine currentLine;
    currentLine.PnId = PnId;
    string tmp;
    if (useFirstLine)
    {
        lexicalCast(currentLine.numOfRepeats, firstLine[0]);
        lexicalCast(currentLine.ratioBf,firstLine[1]);
        lexicalCast(currentLine.ratioAf,firstLine[2]);
        lexicalCast(currentLine.locationShift,firstLine[3]);
        lexicalCast(currentLine.mateEditDist,firstLine[4]);
        lexicalCast(currentLine.purity,firstLine[5]);
        lexicalCast(currentLine.ratioOver20In,firstLine[6]);
    }
    else
    {
        if (is_gz)
        {
            attributeFile_gz >> currentLine.numOfRepeats;
            attributeFile_gz >> currentLine.ratioBf;
            attributeFile_gz >> currentLine.ratioAf;
            attributeFile_gz >> currentLine.locationShift;
            attributeFile_gz >> currentLine.mateEditDist;
            attributeFile_gz >> currentLine.purity;
            attributeFile_gz >> currentLine.ratioOver20In;
            attributeFile_gz >> tmp; //repeatSeqFromRead
        }
        else
        {
            attributeFile >> currentLine.numOfRepeats;
            attributeFile >> currentLine.ratioBf;
            attributeFile >> currentLine.ratioAf;
            attributeFile >> currentLine.locationShift;
            attributeFile >> currentLine.mateEditDist;
            attributeFile >> currentLine.purity;
            attributeFile >> currentLine.ratioOver20In;
            attributeFile >> tmp; //repeatSeqFromRead
        }
    }
    if (!useModelAndLabels)
        currentLine.pValue = 0.95;
    else
        currentLine.pValue = getPval(marker, currentLine);
    if (enoughReads)
        PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].pValueSum += currentLine.pValue;
    //Check if the read is a result of a full motif slippage
    float diff1 = fabs(winner-currentLine.numOfRepeats), diff2 = fabs(second-currentLine.numOfRepeats);
    if (std::min(diff1,diff2)>=0.9)
    {
        PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].fullMotifSlippageSum += currentLine.pValue;
        //Check whether slippage removed or added repeat
        if (diff1 < diff2)
        {
            if (winner > currentLine.numOfRepeats)
                markerToStats[marker].negSlippSum += currentLine.pValue;
            else
                markerToStats[marker].posSlippSum += currentLine.pValue;
        }
        else
        {
            if (second > currentLine.numOfRepeats)
                markerToStats[marker].negSlippSum += currentLine.pValue;
            else
                markerToStats[marker].posSlippSum += currentLine.pValue;
        }
    }
    //Determining the initial label of the read
    if (diff1 <= 0.05 || diff2 <= 0.05)
    {
        currentLine.label = 1;
        ++markerToStats[marker].slippCount.p1.i1;
        if (enoughReads)
            markerToStats[marker].slippCount.p1.i2 += currentLine.pValue;
    }
    else
    {
        if (fabs(currentLine.numOfRepeats - (winner - 1)) <= 0.05 || fabs(currentLine.numOfRepeats - (second - 1)) <= 0.05)
        {
            currentLine.label = 2;
            ++markerToStats[marker].slippCount.p2.i1;
            if (enoughReads)
                markerToStats[marker].slippCount.p2.i2 += currentLine.pValue;
            markerToStats[marker].stepSum += (float)marker.motif.size();
        }
        else
        {
            currentLine.label = -1;
            ++markerToStats[marker].slippCount.p3.i1;
            if (enoughReads)
                markerToStats[marker].slippCount.p3.i2 += currentLine.pValue;
            if (fmod(diff1,1.0)<0.05)
            {
                if (fmod(diff2,1.0)<0.05)
                    markerToStats[marker].stepSum += round(std::min(diff1,diff2) * (float)marker.motif.size());
                else
                    markerToStats[marker].stepSum += round(diff1 * (float)marker.motif.size());
            }
            else
            {
                if (fmod(diff2,1.0)<0.05)
                    markerToStats[marker].stepSum += round(diff2 * (float)marker.motif.size());
                else
                    markerToStats[marker].stepSum += round(std::min(diff1,diff2) * (float)marker.motif.size());
            }
        }
    }
    return currentLine;
}

MakeGenotypesRet makeGenotypes(std::set<float> alleles, float refAllele)
{
    String<Pair<float> > genotypes;
    String<float> alleleString;
    std::set<Pair<float> > genotypeSet;
    std::set<float>::reverse_iterator allelesBegin = alleles.rend();
    float closestToRef = *alleles.rbegin();
    float minDistToRef = fabs(closestToRef-refAllele);
    unsigned refIdx = 0;
    for (std::set<float>::reverse_iterator alleleIt = alleles.rbegin(); alleleIt!=allelesBegin; ++alleleIt)
    {
        appendValue(alleleString, *alleleIt);
        if (fabs(*alleleIt-refAllele)<minDistToRef)
        {
            closestToRef = *alleleIt;
            minDistToRef = fabs(*alleleIt-refAllele);
            refIdx = length(alleleString)-1;
        }
    }
    //Put reference allele at end of alleleString
    erase(alleleString, refIdx);
    appendValue(alleleString, closestToRef);
    for (unsigned i=0; i<length(alleleString); ++i)
    {
        appendValue(genotypes,Pair<float>(alleleString[i],alleleString[i]));
        genotypeSet.insert(Pair<float>(alleleString[i],alleleString[i]));
        if (i == (length(alleleString)-1))
            break;
        for (unsigned j=i+1; j<length(alleleString); ++j)
        {
            appendValue(genotypes,Pair<float>(alleleString[j],alleleString[i]));
            genotypeSet.insert(Pair<float>(alleleString[j],alleleString[i]));
        }
    }
    reverse(genotypes);
    MakeGenotypesRet returnValue;
    returnValue.genotypes = genotypes;
    returnValue.genotypesSet = genotypeSet;
    return returnValue;
}

int findMinIndex(String<long double> probs)
{
    int minIndex = 0;
    long double minValue = probs[0];
    for (unsigned i = 1; i<length(probs); ++i)
    {
        if (probs[i]<minValue)
        {
            minIndex=i;
            minValue=probs[i];
        }
    }
    return minIndex;
}

float dgeom(int diff, double psucc)
{
    if (diff < 0)
        return 0;
    double p = psucc;
    for (int i = 0; i < diff; i++)
    {
        p = p*(1-psucc);
    }
    return p;
}

float dpois(int step, float lambda) {
  if (step < 0)
      return 0;
  float p = exp(-1*lambda);
  for (int i = 0; i < step; i++) {
    p = p*lambda;
    p = p/(i+1);
  }
  return p;
}

Pair<GenotypeInfo, bool> determineGenotype(String<AttributeLine> reads, double markerSlippage, String<Pair<float> > genotypes, int numberOfAlleles, int motifLength, double psucc, double posSlippProb, double negSlippProb)
{
    GenotypeInfo returnValue;
    returnValue.pValueSum = 0;
    returnValue.genotypes = genotypes;
    returnValue.numOfReads = length(reads);
    Pair<float> genotypeToCheck;
    AttributeLine readToCheck;
    String<long double> probs;
    double errorProbSum = 0;
    resize(probs, length(genotypes));
    bool isHomo, enoughDistance = true;
    float diff, diff2, lambda = std::max((double)0.001,markerSlippage), posNegSlipp = 1, posNegSlipp2 = 1;
    int indexOfWinner, indexOfSecond;
    for (unsigned i=0; i<length(genotypes); ++i)
    {
        probs[i] = 0;
        genotypeToCheck = genotypes[i];
        isHomo = genotypeToCheck.i1 == genotypeToCheck.i2;
        for (unsigned j=0; j<length(reads); ++j)
        {
            posNegSlipp = 1;
            posNegSlipp2 = 1;
            readToCheck = reads[j];
            if (i == 0)
            {
                ++returnValue.alleleToFreq[readToCheck.numOfRepeats];
                if (length(reads) >= 10)
                    returnValue.pValueSum += readToCheck.pValue;
            }
            if (isHomo)
            {
                if (readToCheck.numOfRepeats - genotypeToCheck.i1 < -0.9)
                    posNegSlipp = negSlippProb;
                if (readToCheck.numOfRepeats - genotypeToCheck.i1 > 0.9)
                    posNegSlipp = posSlippProb;
                diff = fabs(readToCheck.numOfRepeats - genotypeToCheck.i1);
                if (useGeom)
                    probs[i] += -(double)10*log10(readToCheck.pValue * dgeom(static_cast<int>(roundf((diff-(float)floor(diff))*motifLength)), psucc) * dpois(floor(diff), lambda) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                else
                    probs[i] += -(double)10*log10(readToCheck.pValue * dpois(ceil(diff), lambda) * posNegSlipp + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
            }
            else
            {
                if (readToCheck.numOfRepeats - genotypeToCheck.i1 < -0.9)
                    posNegSlipp = negSlippProb;
                if (readToCheck.numOfRepeats - genotypeToCheck.i1 > 0.9)
                    posNegSlipp = posSlippProb;
                if (readToCheck.numOfRepeats - genotypeToCheck.i2 < -0.9)
                    posNegSlipp2 = negSlippProb;
                if (readToCheck.numOfRepeats - genotypeToCheck.i2 > 0.9)
                    posNegSlipp2 = posSlippProb;
                diff = fabs(readToCheck.numOfRepeats - genotypeToCheck.i1);
                diff2 = fabs(readToCheck.numOfRepeats - genotypeToCheck.i2);
                if (useGeom)
                    probs[i] += -(double)10*log10(readToCheck.pValue * 0.5 * (dgeom(static_cast<int>(roundf((diff-(float)floor(diff))*motifLength)), psucc) * dpois(floor(diff), lambda) * posNegSlipp + dgeom(static_cast<int>(roundf((diff2-(float)floor(diff2))*motifLength)), psucc) * dpois(floor(diff2), lambda) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
                else
                    probs[i] += -(double)10*log10(readToCheck.pValue * 0.5 * (dpois(ceil(diff), lambda) * posNegSlipp + dpois(ceil(diff2), lambda) * posNegSlipp2) + ((double)(1.0-readToCheck.pValue)/(double)numberOfAlleles));
            }
        }
    }
    returnValue.pValues = probs;
    indexOfWinner = findMinIndex(probs);
    if (length(probs)>1)
    {
        String<long double> probsCopy = probs;
        erase(probsCopy,indexOfWinner);
        indexOfSecond = findMinIndex(probsCopy);
        enoughDistance = round(probsCopy[indexOfSecond]-probs[indexOfWinner]) > 10;
    }
    returnValue.genotype = genotypes[indexOfWinner];
    returnValue.pValue = probs[indexOfWinner];
    for (unsigned j=0; j<length(reads); ++j)
    {
       readToCheck = reads[j];
       if ((readToCheck.numOfRepeats != returnValue.genotype.i1) && (readToCheck.numOfRepeats != returnValue.genotype.i2))
        errorProbSum += readToCheck.pValue;
    }
    return Pair<GenotypeInfo, bool>(returnValue, enoughDistance);
}

unsigned getChrLength(string chrom)
{
	if (chrom.compare("chr1") == 0)
		return 248956422;
	if (chrom.compare("chr2") == 0)
		return 242193529;
	if (chrom.compare("chr3") == 0)
		return 198295559;
	if (chrom.compare("chr4") == 0)
		return 190214555;
	if (chrom.compare("chr5") == 0)
		return 181538259;
	if (chrom.compare("chr6") == 0)
		return 170805979;
	if (chrom.compare("chr7") == 0)
		return 159345973;
	if (chrom.compare("chr8") == 0)
		return 145138636;
	if (chrom.compare("chr9") == 0)
		return 138394717;
	if (chrom.compare("chr10") == 0)
		return 133797422;
	if (chrom.compare("chr11") == 0)
		return 135086622;
	if (chrom.compare("chr12") == 0)
		return 133275309;
	if (chrom.compare("chr13") == 0)
		return 114364328;
	if (chrom.compare("chr14") == 0)
		return 107043718;
	if (chrom.compare("chr15") == 0)
		return 101991189;
	if (chrom.compare("chr16") == 0)
		return 90338345;
	if (chrom.compare("chr17") == 0)
		return 83257441;
	if (chrom.compare("chr18") == 0)
		return 80373285;
	if (chrom.compare("chr19") == 0)
		return 58617616;
	if (chrom.compare("chr20") == 0)
		return 64444167;
	if (chrom.compare("chr21") == 0)
		return 46709983;
	if (chrom.compare("chr22") == 0)
		return 50818468;
	else
		return 0;
}

//Write all sorts of info to the header of the vfc file I pass to the function
void makeVcfHeader(VcfFileOut& out, String<string> PnIds, string chrom)
{
    appendValue(contigNames(context(out)), chrom);
    //Add IDs of all PNs to the header
    for (unsigned i = 0; i<length(PnIds); ++i)
        appendValue(sampleNames(context(out)), PnIds[i]);
    unsigned chromLength = getChrLength(chrom);
    string contigString = "<ID=" + chrom + ",length=" + to_string((long long unsigned int)chromLength) + ">";
    //Complicated way of getting todays date
    time_t rawtime;
    tm* timeinfo;
    char buffer [80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,10,"%Y%m%d",timeinfo);
    string date = buffer;
    VcfHeader header;
    appendValue(header, VcfHeaderRecord("fileformat", "VCFv4.2"));
    appendValue(header, VcfHeaderRecord("fileDate", date));
    appendValue(header, VcfHeaderRecord("source", "PopSTR"));
    appendValue(header, VcfHeaderRecord("source_bin", "/odinn/tmp/snaedisk/190423_popSTRfreeze/binary/computeReadAttributes"));
    appendValue(header, VcfHeaderRecord("source_bin", "/odinn/tmp/snaedisk/190423_popSTRfreeze/binary/computePnSlippage"));
    appendValue(header, VcfHeaderRecord("source_bin", "/odinn/tmp/snaedisk/190423_popSTRfreeze/binary/msGenotyper"));
    appendValue(header, VcfHeaderRecord("reference", "/odinn/data/reference/Homo_sapiens-deCODE-hg38/Sequence/WholeGenomeFasta/genome.fa"));
    appendValue(header, VcfHeaderRecord("contig", contigString));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=RefLen,Number=A,Type=Integer,Description=\"Length of the reference allele\">"));
    appendValue(header, VcfHeaderRecord("INFO", "<ID=Motif,Number=1,Type=String,Description=\"Microsatellite repeat motif\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">"));
    appendValue(header, VcfHeaderRecord("FORMAT", "<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods\">"));
    writeHeader(out, header);
}

//Clear a stringstream and a string, I use this a lot in fillRecordMarker/Pn, so a function seemed appropriate
void stringClear(stringstream& ss, CharString& str)
{
    ss.str("");
    ss.clear();
    clear(str);
}

Pair<int> findGenotypeIndices(GenotypeInfo& genotype, std::set<float>& allelesAtMarker)
{
    Pair<int> gtIdxs = Pair<int>(0,0);
    unsigned idx = 1;
    float currentAllele;
    std::set<float>::iterator allEnd = allelesAtMarker.end();
    for (std::set<float>::iterator allIt = allelesAtMarker.begin(); allIt!=allEnd; ++allIt)
    {
        currentAllele = *allIt;
        if (currentAllele == genotype.genotype.i1)
            gtIdxs.i1 = idx;
        if (currentAllele == genotype.genotype.i2)
            gtIdxs.i2 = idx;
        idx++;
    }
    return gtIdxs;
}

VcfRecord fillRecordMarker(Marker marker, std::set<float> allelesAtThisMarker)
{
    int refLength = marker.end - (marker.start-1);
    VcfRecord record;
    record.rID = 0;
    record.beginPos = marker.start-1;
    stringstream ss;
    ss << record.beginPos+1;
    CharString str = ss.str();
    record.id = marker.chrom + ":" + toCString(str) + ":M";
    //Get ref seq using crop_fasta and set as record.ref
    string refSeq;
    FILE * stream;
    char buffer[256];
    stringstream cmdStream;
    cmdStream << "/odinn/users/snaedisk/bin/crop_fasta /odinn/data/reference/Homo_sapiens-deCODE-hg38/Sequence/WholeGenomeFasta/genome.fa " << marker.chrom << ":" << toCString(str) << "-" << marker.end;
    string cmd = cmdStream.str();
    cmd.append(" 2>&1");
    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
        if (fgets(buffer, 256, stream) != NULL)
        {
            buffer[strcspn(buffer, "\n")] = 0;
            refSeq.append(buffer);
        }
        pclose(stream);
    }
    //Set reference allele and delete it from allele set
    stringClear(ss,str);
    record.ref = refSeq;
    allelesAtThisMarker.erase(marker.refRepeatNum);
    record.qual = record.MISSING_QUAL();
    record.info = "RefLen=";
    ss << refLength;
	str = ss.str();
	append(record.info,str);
	stringClear(ss,str);
    append(record.info,";");
    append(record.info,"Motif=");
    append(record.info,marker.motif);
    record.format = "GT:AD:DP:PL";
    //Loop over all alternative alleles and add them to ALT field in record
    float currentAllele;
    std::set<float>::iterator allEnd = allelesAtThisMarker.end();
    for (std::set<float>::iterator allIt = allelesAtThisMarker.begin(); allIt!=allEnd; ++allIt)
    {
        append(record.alt, "<");
        currentAllele = *allIt;
        ss << currentAllele;
        str = ss.str();
        append(record.alt, str);
        append(record.alt, ">,");
        stringClear(ss,str);
    }
    if (allelesAtThisMarker.size() > 0)
        eraseBack(record.alt);
    else
        record.alt = ".";
    return record;
}

//Fill up the PN-specific fields in the VCF-record.
void fillRecordPn(GenotypeInfo genotype, VcfRecord& record, MakeGenotypesRet genotypesAtThisMarker, std::set<float>& allelesAtMarker, Marker& marker)
{
    unsigned indexOfCurrGt;
    while (genotypesAtThisMarker.genotypesSet.find(genotype.genotype)== genotypesAtThisMarker.genotypesSet.end() && length(genotype.genotypes)>1)
    {
        indexOfCurrGt = findMinIndex(genotype.pValues);
        erase(genotype.pValues,indexOfCurrGt);
        erase(genotype.genotypes,indexOfCurrGt);
        indexOfCurrGt = findMinIndex(genotype.pValues);
        genotype.genotype = genotype.genotypes[indexOfCurrGt];
        genotype.pValue = genotype.pValues[indexOfCurrGt];
    }
    float refAllele = marker.refRepeatNum;
    //Remove ref allele from set.
    allelesAtMarker.erase(refAllele);
    stringstream ss;
    CharString str = ss.str();
    CharString gtInfo; //First I make the string containing the genotype info
    Pair<int> gtIdxs = findGenotypeIndices(genotype, allelesAtMarker);
    ss << gtIdxs.i1;
    str = ss.str();
    append(gtInfo,str);
    stringClear(ss,str);
    append(gtInfo,"/");
    ss << gtIdxs.i2;
    str = ss.str();
    append(gtInfo,str);
    stringClear(ss,str);
    append(gtInfo,":");
    //Add alleleDepth count for the ref allele
    if (genotype.alleleToFreq.count(refAllele)==0)
        append(gtInfo,"0,"); //No reads supporting ref allele
    else
    {
        ss << genotype.alleleToFreq[refAllele];
        str = ss.str();
        append(gtInfo,str);
        stringClear(ss,str);
        append(gtInfo,",");
    }
    //Remove ref allele from set before looping over
    allelesAtMarker.erase(refAllele);
    //Add allele depth counts for other alleles
    for (auto const & allele : allelesAtMarker)
    {
        if (genotype.alleleToFreq.count(allele) == 0)
            append(gtInfo,"0,");
        else
        {
            ss << genotype.alleleToFreq[allele];
            str = ss.str();
            append(gtInfo,str);
            stringClear(ss,str);
            append(gtInfo,",");
        }
    }
    //put ref allele back in set after looping over
    allelesAtMarker.insert(refAllele);
    unsigned possibleGenotypes = 0;
    for (unsigned i=0; i<allelesAtMarker.size(); ++i)
    {
        possibleGenotypes += i+1;
    }
    if (length(genotypesAtThisMarker.genotypes) != possibleGenotypes)
        cerr << "ATTENTION: number of genotypes considered is wrong. It should be: " << possibleGenotypes << " but it is:" << length(genotypesAtThisMarker.genotypes) << endl;
    eraseBack(gtInfo);
    append(gtInfo,":");
    //Add readDepth to gtInfo
    ss << genotype.numOfReads;
    str = ss.str();
    append(gtInfo,str);
    stringClear(ss,str);
    append(gtInfo,":");
    Pair<float> genotypeToLookFor;
    Pair<float> genotypeToCompare;
    int index;
    long double numerator;
    long double denominator;
    int pl;
    for (unsigned i=0; i<length(genotypesAtThisMarker.genotypes); ++i)
    {
        index = -1;
        genotypeToLookFor = genotypesAtThisMarker.genotypes[i];
        if (genotypeToLookFor == genotype.genotype)
        {
            ss << 0;
            str = ss.str();
            append(gtInfo,str);
        }
        else
        {
            for (unsigned j=0; j<length(genotype.genotypes); ++j)
            {
                genotypeToCompare=genotype.genotypes[j];
                if (genotypeToCompare == genotypeToLookFor)
                {
                    index = j;
                    break;
                }
            }
            if (index == -1)
            {
                ss << 255;
                str = ss.str();
                append(gtInfo,str);
            }
            else
            {
                numerator = genotype.pValues[index];
                denominator = genotype.pValue;
                pl = round(numerator-denominator);
                pl = std::min(255,pl);
                if (pl<0)
                    pl = 255;
                ss << pl;
                str = ss.str();
                append(gtInfo,str);
            }
        }
        stringClear(ss,str);
        append(gtInfo,",");
    }
    eraseBack(gtInfo);
    //When the genotype info string is ready I append it to the stringset of charstrings
    appendValue(record.genotypeInfos, gtInfo);
}

double computeAlleleDist(String<Pair<float> > genotypes, map<float,int> allelesToFreq, int pnsAtMarker)
{
    float genotype1, genotype2;
    double freq1, freq2;
    double currVal;
    double totalSum = 0;
    int indicator;
    for (unsigned i=0; i<length(genotypes); ++i)
    {
        indicator = 0;
        genotype1 = genotypes[i].i1;
        genotype2 = genotypes[i].i2;
        float n1 = fabs(genotype1-genotype2);
        float n2 = ceil(fabs(genotype1-genotype2));
        float n3 = n1-n2;
        if (n3 == 0)
            indicator = 1;
        freq1 = (double)allelesToFreq[genotype1]/(double)(2*pnsAtMarker);
        freq2 = (double)allelesToFreq[genotype2]/(double)(2*pnsAtMarker);
        currVal = (double)indicator*freq1*freq2;
        if (genotype1 != genotype2)
            currVal = 2*currVal;
        totalSum += currVal;
    }
    return totalSum;
}

//Count number of words in a sentence, use to parse input from attribute file
Pair<int, String<string> > countNumberOfWords(string sentence)
{
    int numberOfWords = 0;
    String<string> words;
    resize(words, 9);
    int currentWordLength;
    if (!isspace(sentence[0]))
    {
        numberOfWords++;
        words[0] = sentence[0];
    }

    for (unsigned i = 1; i < sentence.length(); i++)
    {
        if ((!isspace(sentence[i])) && (isspace(sentence[i-1])))
        {
            numberOfWords++;
            words[numberOfWords-1] = sentence[i];
        }
        else
        {
            if (!isspace(sentence[i]))
                words[numberOfWords-1].push_back(sentence[i]);
        }
    }

    resize(words, numberOfWords);
    return Pair<int, String<string> >(numberOfWords, words);
}

String<Pair<int, string> > readMarkers(CharString & markerFile)
{
    ifstream markerStream(toCString(markerFile));
    if(markerStream.fail())
        cout << "Unable to locate markerList @ " << markerFile << endl;
    String<Pair<int, string> > markerList;
    string motif;
    int startPos;
    while (true)
    {
        markerStream >> startPos;
        markerStream >> motif;
        if (markerStream.eof() || motif.length() == 0)
            break;
        Pair<int, string> currPair = Pair<int, string>(startPos, motif);
        appendValue(markerList, currPair);
    }
    markerStream.close();
    cout << "Finished reading " << length(markerList) << " markers.\n";
    return markerList;
}

//Read in data from markerSlippageFile
void readMarkerSlippage(CharString markerSlippageFile, CharString iterationNumber, CharString regressionModelDirectory, CharString intervalIndex)
{
    string prevItNum;
    int itNum = 0;
    lexicalCast(itNum, iterationNumber);
    --itNum;
    prevItNum = to_string(itNum);
    append(markerSlippageFile, prevItNum);
    append(markerSlippageFile, "_");
    append(markerSlippageFile, intervalIndex);
    ifstream markerSlippageIn(toCString(markerSlippageFile));
    double currMarkSlipp, currMarkStutt, currPosSlippProb, currNegSlippProb;
    int nPns, nAlleles;
    Marker currMarker;
    CharString currMarkerModelDir = regressionModelDirectory;
    while (true)
    {
        markerSlippageIn >> currMarker.chrom;
        markerSlippageIn >> currMarker.start;
        markerSlippageIn >> currMarker.end;
        markerSlippageIn >> currMarker.motif;
        markerSlippageIn >> currMarker.refRepeatNum;
        markerSlippageIn >> currMarkSlipp;
        markerSlippageIn >> nPns;
        markerSlippageIn >> nAlleles;
        markerSlippageIn >> currMarkStutt;
        markerSlippageIn >> currPosSlippProb;
        markerSlippageIn >> currNegSlippProb;
        if (markerSlippageIn.eof())
            break;
        markerToStats[currMarker].nAlleles = nAlleles;
        markerToStats[currMarker].stutter = currMarkStutt;
        markerToStats[currMarker].slippage = currMarkSlipp;
        append(currMarkerModelDir, "/model_");
        append(currMarkerModelDir, to_string(currMarker.start));
        append(currMarkerModelDir, "_");
        append(currMarkerModelDir, currMarker.motif);
        const char *model_in_file = toCString(currMarkerModelDir);
        markerToStats[currMarker].regressionModel = load_model(model_in_file);
        currMarkerModelDir = regressionModelDirectory;
    }
    cout << "Finished reading marker slippage." << endl;
    markerSlippageIn.close();
}

map<string, Pair<double, int> > readPnSlippage(CharString pnSlippagePath, CharString iterationNumber, String<string> & PnIds)
{
    string PnId;
    int nMarkers;
    double currPnSlipp;
    string prevItNum;
    int itNum = 0;
    lexicalCast(itNum, iterationNumber);
    --itNum;
    prevItNum = to_string(itNum);
    map<string, Pair<double, int> > pnToSize;
    append(pnSlippagePath, "_");
    append(pnSlippagePath, prevItNum);
    ifstream pnSlippageFile(toCString(pnSlippagePath));
    while (true)
    {
        pnSlippageFile >> PnId;
        pnSlippageFile >> currPnSlipp;
        pnSlippageFile >> nMarkers;
        if (pnSlippageFile.eof() || PnId.length() == 0)
            break;
        pnToSize[PnId].i1= currPnSlipp;
        pnToSize[PnId].i2 = nMarkers;
        appendValue(PnIds, PnId);
    }
    cout << "Finished reading pn Slippage, number of pns:" << pnToSize.size() << endl;
    assert (pnToSize.size() == length(PnIds));
    pnSlippageFile.close();
    return pnToSize;
}

inline bool exists(const std::string& name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

Pair<double, int> estimateSlippage(String<string> PnIds, map<Pair<string,Marker>, GenotypeInfo> PnAndMarkerToGenotype, Marker marker, CharString currItNum, map<string, Pair<double, int> > pnToSize)
{
    int nMissing = 0, nAvailable;
    vector<double> weights;
    vector<double> slippFragments;
    double currPnSlipp, currPvalSum = 0, weightSum = 0, fullMotifSlippageSum = 0, result, currMarkSlipp;
    if (atoi(toCString(currItNum)) == 1)
    {
        for (unsigned i = 0; i<length(PnIds); ++i)
        {
            if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(PnIds[i], marker)) == 0) || (PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].pValueSum == 0))
                continue;
            fullMotifSlippageSum += PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].fullMotifSlippageSum;
            currPvalSum += PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i],marker)].pValueSum;
        }
        currMarkSlipp = 0.5*(fullMotifSlippageSum/currPvalSum);
    }
    else
        currMarkSlipp = markerToStats[marker].slippage;
    for (unsigned i = 0; i<length(PnIds); ++i)
    {
        if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(PnIds[i], marker)) == 0) || (PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].pValueSum == 0))
        {
            weights.push_back(0);
            ++nMissing;
            continue;
        }
        if (pnToSize[PnIds[i]].i1 == 0)
            currPnSlipp = 0.001;
        else
            currPnSlipp = pnToSize[PnIds[i]].i1;
        currPvalSum = PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i],marker)].pValueSum;
        weights.push_back(currPvalSum/((currPnSlipp+currMarkSlipp)*(1-(currPnSlipp+currMarkSlipp))));
    }
    weightSum = accumulate(weights.begin(),weights.end(),0.0);
    nAvailable = length(PnIds) - nMissing;
    for (unsigned i = 0; i<length(PnIds); ++i)
    {
        if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(PnIds[i], marker)) == 0) || (PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].pValueSum == 0))
        {
            slippFragments.push_back(0);
            continue;
        }
        if (pnToSize[PnIds[i]].i1 == 0)
            currPnSlipp = 0.001;
        else
            currPnSlipp = pnToSize[PnIds[i]].i1;
        fullMotifSlippageSum = PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i], marker)].fullMotifSlippageSum;
        currPvalSum = PnAndMarkerToGenotype[Pair<string,Marker>(PnIds[i],marker)].pValueSum;
        slippFragments.push_back((weights[i]/weightSum)*((fullMotifSlippageSum/currPvalSum)-currPnSlipp));
    }
    result = std::max(0.0,accumulate(slippFragments.begin(),slippFragments.end(),0.0));
    cout << "Slippage rate: " << result << endl;
    return Pair<double, int>(result, nAvailable);
}

Pair<bool> genotypeIsConfident(GenotypeInfo& genotypeInfo)
{
    bool A1ok = false, A2ok = false;
    if (genotypeInfo.alleleToFreq[genotypeInfo.genotype.i1]>=3 && (double)genotypeInfo.alleleToFreq[genotypeInfo.genotype.i1]/(double)genotypeInfo.numOfReads >=0.25)
        A1ok = true;
    if (genotypeInfo.alleleToFreq[genotypeInfo.genotype.i2]>=3 && (double)genotypeInfo.alleleToFreq[genotypeInfo.genotype.i2]/(double)genotypeInfo.numOfReads >=0.25)
        A2ok = true;
    return Pair<bool>(A1ok, A2ok);
}

//For storing command line arguments
struct MsGenotyperOptions
{
    CharString attDirChromNum, pnSlippageFile, markerList, intervalIndex, markerSlippageFile, regressionModelDirectory, iterationNumber, vcfOutputDirectory, vcfFileName;
    unsigned firstPnIdx;
} ;

ArgumentParser::ParseResult parseCommandLine(MsGenotyperOptions & options, int argc, char const ** argv)
{
    ArgumentParser parser("msGenotyper");
    setShortDescription(parser, "Microsatellite genotyper");
    setVersion(parser, "1.3");
    setDate(parser, "June 2019");
    addUsageLine(parser, "\\fI-AD\\fP attributesDirectory/chromNum \\fI-PNS\\fP pnSlippageFile \\fI-ML\\fP markerList \\fI-I\\fP intervalIndex \\fI-MS\\fP markerSlippageFile \\fI-MD\\fP regressionModelDirectory \\fI-IN\\fP iterationNumber \\fI-FP\\fP firstPnIdx \\fI-VD\\fP vcfOutputDirectory \\fI-VN\\fP vcfFileName  ");
    addDescription(parser, "Performs genptyping for all PNs in the pnSlippageFile over all markers in the markerFile. The genotypes are written to a file specified by the user.");

    addOption(parser, ArgParseOption("AD", "attributesDirectory/chromNum", "Path to attributes files for the markers being genotyped.", ArgParseArgument::INPUT_FILE, "IN-DIR"));
    setRequired(parser, "attributesDirectory/chromNum");

    addOption(parser, ArgParseOption("PNS", "pnSlippageFile", "A file containing slippage rates for the pns to be genotyped.", ArgParseArgument::INPUT_FILE, "IN-FILE"));
    setRequired(parser, "pnSlippageFile");

    addOption(parser, ArgParseOption("ML", "markerList", "List of markers to genotype.", ArgParseArgument::INPUT_FILE, "IN-FILE"));
    setRequired(parser, "markerList");

    addOption(parser, ArgParseOption("I", "intervalIndex", "Index of the interval being processed", ArgParseArgument::STRING, "INDEX"));
    setRequired(parser, "intervalIndex");

    addOption(parser, ArgParseOption("MS", "markerSlippageFile", "A file containing estimated slippage rates for the microsatellites.", ArgParseArgument::OUTPUT_FILE, "OUT-FILE"));
    setRequired(parser, "markerSlippageFile");

    addOption(parser, ArgParseOption("MD", "regressionModelDirectory", "A directory where logistic regression models for all markers in the markerList are stored", ArgParseArgument::INPUTPREFIX, "IN-DIR"));
    setRequired(parser, "regressionModelDirectory");

    addOption(parser, ArgParseOption("IN", "iterationNumber", "Index of the iteration being performed, 1-based.", ArgParseArgument::STRING, "INDEX"));
    setRequired(parser, "iterationNumber");

    addOption(parser, ArgParseOption("FP", "firstPnIdx", "Index of first Pn in pnList within the attributeFile.", ArgParseArgument::INTEGER, "INTEGER"));
    setRequired(parser, "firstPnIdx");

    //Used if this is the final iteration
    addOption(parser, ArgParseOption("VD", "vcfOutputDirectory", "A directory to write the vcf file to.", ArgParseArgument::OUTPUT_FILE, "OUT-DIR"));

    addOption(parser, ArgParseOption("VN", "vcfFileName", "Name of vcf output file.", ArgParseArgument::OUTPUT_FILE, "OUT-FILE"));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    if (res != ArgumentParser::PARSE_OK)
        return res;

    getOptionValue(options.attDirChromNum, parser, "attributesDirectory/chromNum");
    append(options.attDirChromNum, "/");
    getOptionValue(options.pnSlippageFile, parser, "pnSlippageFile");
    getOptionValue(options.markerList, parser, "markerList");
    getOptionValue(options.intervalIndex, parser, "intervalIndex");
    getOptionValue(options.markerSlippageFile, parser, "markerSlippageFile");
    getOptionValue(options.regressionModelDirectory, parser, "regressionModelDirectory");
    getOptionValue(options.iterationNumber, parser, "iterationNumber");
    getOptionValue(options.firstPnIdx, parser, "firstPnIdx");

    if (isSet(parser,"vcfOutputDirectory"))
    {
        getOptionValue(options.vcfOutputDirectory, parser, "vcfOutputDirectory");
        getOptionValue(options.vcfFileName, parser, "vcfFileName");
    }
    else
    {
        options.vcfOutputDirectory = ".";
        options.vcfFileName = "deleteMe";
    }

    return ArgumentParser::PARSE_OK;
}

map<string, Pair<float> > readLabels(CharString modelAndLabelDir, Pair<int, string> marker, CharString iterationNumber)
{
    map<string, Pair<float> > pnToLabels;
    append(modelAndLabelDir, "/");
    append(modelAndLabelDir, to_string(marker.i1));
    append(modelAndLabelDir, "_");
    append(modelAndLabelDir, marker.i2);
    string prevItNum;
    int itNum = 0;
    lexicalCast(itNum, iterationNumber);
    --itNum;
    prevItNum = to_string(itNum);
    append(modelAndLabelDir, prevItNum);
    ifstream labelFile(toCString(modelAndLabelDir));
    while (!labelFile.eof())
    {
        string PN_ID;
        float A1, A2;
        labelFile >> PN_ID;
        if (PN_ID.length() == 0)
            break;
        labelFile >> A1 >> A2;
        Pair<float> labels = Pair<float>(A1, A2);
        pnToLabels[PN_ID] = labels;
    }
    labelFile.close();
    return pnToLabels;
}

long int readOffSets(ifstream & attsFile, unsigned firstPnIdx)
{
    long int offset = 0;
    for (unsigned i = 1; i<=firstPnIdx; ++i)
    {
        attsFile >> offset;
        if (offset == -69)
            return 0;
    }
    while (offset == 0)
        attsFile >> offset;
    return offset;
}

int main(int argc, char const ** argv)
{
    //Check arguments.
    MsGenotyperOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    //Read marker and pn data
    String<string> PnIds;
    map<string, Pair<double, int> > pnToSize = readPnSlippage(options.pnSlippageFile, options.iterationNumber, PnIds); //pnslippage
    String<Pair<int, string> > markersForGenotyping = readMarkers(options.markerList); //list of markers to genotype
    bool loadModAndLab = true;
    if (atoi(toCString(options.iterationNumber)) > 1) //previous estimates and regression models, if available
        readMarkerSlippage(options.markerSlippageFile, options.iterationNumber, options.regressionModelDirectory, options.intervalIndex);
    else
        loadModAndLab = false;

    //Variable declarations and initializations
    string PnId, chrom, motif, nextWord, finalItNum = "5";
    int start, end, numberOfReads;
    bool writeVcf = false, enoughReads = true;
    if (options.iterationNumber == finalItNum)
        writeVcf = true;
    float refRepeatNum, winner, second;

    //Make stream to write markerSlippage estimations to
    append(options.markerSlippageFile, options.iterationNumber);
    append(options.markerSlippageFile, "_");
    append(options.markerSlippageFile, options.intervalIndex);
    ofstream markerSlippageOut(toCString(options.markerSlippageFile));

    //Map from marker to all reads covering it
    map<Marker, String<AttributeLine> > mapPerMarker;
    //Count the total number of alleles in the population for each marker -- by checking the size of the set
    map<Marker, std::set<float> > markerToAlleles;
    //Map to store current genotype(and lots of other things) of a person for each marker maps from pnId and Marker-struct to GenotypeInfo struct
    map<Pair<string,Marker>, GenotypeInfo> PnAndMarkerToGenotype;

    Marker marker;
    string nextLine, temp;
    AttributeLine currentLine;
    Pair<int, String<string> > numberOfWordsAndWords;
    map<string, Pair<float> > pnToLabels;

    //Iterate over all markers in the markerList and read from their attribute files
    int nProcessedMarkers = 0;
    CharString attributePath = options.attDirChromNum;
    for (unsigned i=0; i<length(markersForGenotyping); ++i)
    {
        //Read labelPath for marker, if available
        if (loadModAndLab)
            pnToLabels = readLabels(options.regressionModelDirectory, markersForGenotyping[i], options.iterationNumber);
        //Make path to attributefile for current marker
        append(attributePath, to_string(markersForGenotyping[i].i1));
        append(attributePath, "_");
        append(attributePath, markersForGenotyping[i].i2);
        bool is_gz = false;
        io::filtering_istream attributeFile_gz;
        ifstream attributeFile;
        if (!ifstream(toCString(attributePath)))
        {
            is_gz = true;
            append(attributePath, ".gz");
            std::ifstream attributeFileGz(toCString(attributePath), std::ios_base::in | std::ios_base::binary);
            attributeFile_gz.push(io::gzip_decompressor());
            attributeFile_gz.push(attributeFileGz);
            std::getline (attributeFile_gz, nextLine);
            PnId = nextLine;
        }
        else
        {
            attributeFile.open(toCString(attributePath));
            long int offset = readOffSets(attributeFile, options.firstPnIdx);
            if (offset != 0)
                attributeFile.seekg(offset);
            else
                continue;
        }
        ++nProcessedMarkers;
        if (nProcessedMarkers % 100==0)
            cout << "Working on marker number: " << nProcessedMarkers << endl;
        unsigned pnsFound = 0;
        while (is_gz ? std::getline (attributeFile_gz, nextLine) : std::getline (attributeFile, nextLine))
        {
            if (nextLine.length() == 0)
                continue;
            if (pnsFound >= length(PnIds))
                break;
            numberOfWordsAndWords = countNumberOfWords(nextLine);
            if (numberOfWordsAndWords.i1 == 1)
            {
                //If I have passed the last pn to genotype
                if (nextLine > pnToSize.rbegin()->first)
                    break;
                //Is this PN on my list?
                if (pnToSize.count(nextLine) != 0)
                {
                    ++pnsFound;
                    PnId = nextLine;
                }
                //If it's not I have to read past it.
                else
                {
                    if (is_gz)
                    {
                        attributeFile_gz >> temp;
                        attributeFile_gz >> temp;
                        attributeFile_gz >> temp;
                        attributeFile_gz >> temp;
                        attributeFile_gz >> temp;
                        attributeFile_gz >> numberOfReads;
                        attributeFile_gz >> temp;
                        attributeFile_gz >> temp;
                        for (unsigned i = 0; i <= numberOfReads; ++i)
                            getline (attributeFile_gz,nextLine);
                    }
                    else
                    {
                        attributeFile >> temp;
                        attributeFile >> temp;
                        attributeFile >> temp;
                        attributeFile >> temp;
                        attributeFile >> temp;
                        attributeFile >> numberOfReads;
                        attributeFile >> temp;
                        attributeFile >> temp;
                        for (unsigned i = 0; i <= numberOfReads; ++i)
                            getline (attributeFile,nextLine);
                    }
                }
                continue;
            }
            if (numberOfWordsAndWords.i1 == 8)
            {
                if (numberOfWordsAndWords.i2[0].substr(0,3) == "chr")
                {
                    marker.chrom = numberOfWordsAndWords.i2[0];
                    lexicalCast(marker.start, numberOfWordsAndWords.i2[1]);
                    lexicalCast(marker.end, numberOfWordsAndWords.i2[2]);
                    marker.motif = numberOfWordsAndWords.i2[3];
                    lexicalCast(marker.refRepeatNum, numberOfWordsAndWords.i2[4]);
                    lexicalCast(numberOfReads,numberOfWordsAndWords.i2[5]);
                    if (loadModAndLab)
                    {
                        winner = pnToLabels[PnId].i1;
                        second = pnToLabels[PnId].i2;
                    }
                    else
                    {
                        lexicalCast(winner, numberOfWordsAndWords.i2[6]);
                        lexicalCast(second, numberOfWordsAndWords.i2[7]);
                    }
                    markerToStats[marker].alleles.insert(winner);
                    markerToStats[marker].alleles.insert(second);
                    continue;
                }
                else
                {
                    if (numberOfReads < 10)
                    {
                        PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].pValueSum = 0;
                        PnAndMarkerToGenotype[Pair<string,Marker>(PnId, marker)].fullMotifSlippageSum = 0;
                        enoughReads = false;
                    }
                    for (unsigned i = 0; i < numberOfReads; ++i)
                    {
                        if (i == 0)
                            currentLine = parseNextLine(winner, second, is_gz, attributeFile, attributeFile_gz, marker, PnId, PnAndMarkerToGenotype, numberOfWordsAndWords.i2, true, loadModAndLab, enoughReads);
                        else
                            currentLine = parseNextLine(winner, second, is_gz, attributeFile, attributeFile_gz, marker, PnId, PnAndMarkerToGenotype, numberOfWordsAndWords.i2, false, loadModAndLab, enoughReads);
                        appendValue(mapPerMarker[marker],currentLine);
                    }
                    enoughReads = true;
                    continue;
                }
            }
            if (numberOfWordsAndWords.i1 != 1 && numberOfWordsAndWords.i1 != 8)
                cerr << "Format error in attribute file!" << endl;
        }
        attributePath = options.attDirChromNum;
    }
    pnToLabels.clear();
    chrom = marker.chrom;
    cout << "Reading data from input complete." << endl;

    //Open vcf stream and make header (have to do this every time because otherwise ....segfault)
    VcfFileOut out;
    append(options.vcfOutputDirectory, "/");
    append(options.vcfOutputDirectory, options.vcfFileName);
    append(options.vcfOutputDirectory, "_");
    append(options.vcfOutputDirectory, options.intervalIndex);
    append(options.vcfOutputDirectory, ".vcf");
    ofstream outputFile(toCString(options.vcfOutputDirectory));
    open(out, outputFile);
    makeVcfHeader(out, PnIds, chrom);

    cout << "Number of markers: " << mapPerMarker.size() << endl;

    double *prob_estimates = NULL;
    double predict_label;
    int PnsAtMarker;
    unsigned numOfAlleles;
    //Initialize parameter object for logistic regression
    param.solver_type = 0;
    param.C = 1;
    param.eps = 0.1;
    param.p = 0.1;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;
    //Initialize problem objects for logistic regression training and predicting
    prob.n = 6;
    probBig.n = 6;
    prob.bias = bias;
    probBig.bias = bias;
    //Map to store alleles present in each individual at a given marker. Is cleared for each marker.
    map<string, std::set<float> > PnToAlleles;
    //Store set of all alleles in population for marker being considered and use it to generate all possible genotypes and compute their p-values
    std::set<float> allelesAtMarker;
    std::set<float> allelesToConsider;
    String<Pair<float> > genotypesToConsider;
    float currAllele;
    //Map to store a map from alleles to their frequencies in the population for each marker, used for estimating the probability that the distance between alleles at the marker is an integer and number of PNs at each marker.
    map<Marker, Pair<map<float,int>,int> > markerToAlleleFreqs;
    int z = 1, nAvailable;

    //Stuff for vcf file
    GenotypeInfo genotype;
    string thisPn;
    Marker thisMarker;
    VcfRecord record;
    stringstream ss;
    CharString str;
    MakeGenotypesRet genotypesAtThisMarker;
    double alleleDistance, posSlippProb, negSlippProb;
    Pair<double, int> slippAndNavail;

    //Loop over map from Marker to string<AttributeLine> and train model for each marker and use it to determine genotype
    map<Marker, String<AttributeLine> >::iterator it = mapPerMarker.begin();
    while (mapPerMarker.size()>0)
    {
        thisMarker = it->first;
        PnsAtMarker = 1;
        String<AttributeLine>& currentMarker = it->second;
        probBig.l = length(currentMarker);
        cout << "Starting marker number: " << z << " with start coordinate: " << thisMarker.start << endl;
        double geomP = 1/(fmod(markerToStats[thisMarker].stepSum/(float)length(currentMarker),1.0)+1);
        if (markerToStats[thisMarker].posSlippSum + markerToStats[thisMarker].negSlippSum > 0)
        {
            if (markerToStats[thisMarker].posSlippSum == 0.0)
                markerToStats[thisMarker].posSlippSum = 0.01 * markerToStats[thisMarker].negSlippSum;
            if (markerToStats[thisMarker].negSlippSum == 0.0)
                markerToStats[thisMarker].negSlippSum = 0.01 * markerToStats[thisMarker].posSlippSum;
            posSlippProb = markerToStats[thisMarker].posSlippSum/(markerToStats[thisMarker].posSlippSum + markerToStats[thisMarker].negSlippSum);
            negSlippProb = markerToStats[thisMarker].negSlippSum/(markerToStats[thisMarker].posSlippSum + markerToStats[thisMarker].negSlippSum);
        }
        else
        {
            markerToStats[thisMarker].posSlippSum = 1.0;
            markerToStats[thisMarker].negSlippSum = 1.0;
            posSlippProb = negSlippProb = 0.5;
        }
        allelesAtMarker = markerToStats[thisMarker].alleles;
        numOfAlleles = allelesAtMarker.size();
        markerToStats[thisMarker].nAlleles = numOfAlleles;
        markerToStats[thisMarker].stutter = geomP;
        markerToStats[thisMarker].alleles.clear();
        markerToAlleleFreqs[thisMarker].i1.clear();
        //Estimate marker slippage and update in markerToStats map
        slippAndNavail = estimateSlippage(PnIds, PnAndMarkerToGenotype, thisMarker, options.iterationNumber, pnToSize);
        markerToStats[thisMarker].slippage = std::max((double)0.0,slippAndNavail.i1);
        nAvailable = slippAndNavail.i2;
        //Reads with label 2 are not included in training
        prob.l = length(currentMarker) - markerToStats[thisMarker].slippCount.p2.i1;
        int idx = 0;
        prob.y = Malloc(double,prob.l);
        prob.x = (feature_node **) malloc(prob.l * sizeof(feature_node *));
        probBig.x = (feature_node **) malloc(probBig.l * sizeof(feature_node *));
        PnId = currentMarker[0].PnId;
        for (unsigned i = 0; i<length(currentMarker); ++i)
        {
            currentLine = currentMarker[i];
            if (currentLine.PnId != PnId)
            {
                PnId = currentLine.PnId;
                ++PnsAtMarker;
            }
            PnToAlleles[currentLine.PnId].insert(currentLine.numOfRepeats);
            probBig.x[i] = (feature_node *) malloc(10 * sizeof(feature_node));
            fillProblemX(i,currentLine,probBig);
            //Just include reads with label = 1 or label = -1 in logistic regression training
            if (currentLine.label == 1 || currentLine.label == -1)
            {
                prob.x[idx] = (feature_node *) malloc(10 * sizeof(feature_node));
                fillProblemX(idx,currentLine,prob);
                prob.y[idx] = currentLine.label;
                ++idx;
            }

        }
        const char *error_msg;
        error_msg = check_parameter(&prob,&param);
        if (error_msg != NULL)
            cout << "Error message: " << error_msg << endl;
        //Train logistic regression model
        model_ = train(&prob, &param);
        prob_estimates = (double *) malloc(2*sizeof(double));
        Pair<GenotypeInfo, bool> changed;
        PnId = currentMarker[0].PnId;
        String<AttributeLine> reads;
        for (unsigned i = 0; i < length(currentMarker); ++i)
        {
            currentLine = currentMarker[i];
            if (currentLine.PnId != PnId)
            {
                //Need to reset this set for each PN so I can add alleles from their reads
                allelesToConsider = allelesAtMarker;
                //Add alleles in reads from this PN to the set of alleles present at the marker.
                std::set<float>::iterator pnAllsEnd = PnToAlleles[PnId].end();
                for (std::set<float>::iterator pnAlls = PnToAlleles[PnId].begin(); pnAlls!=pnAllsEnd; ++pnAlls)
                {
                    currAllele = *pnAlls;
                    allelesToConsider.insert(currAllele);
                }
                genotypesToConsider = makeGenotypes(allelesToConsider, thisMarker.refRepeatNum).genotypes;
                //make decision about genotype for PnId at the current marker.
                changed = determineGenotype(reads, markerToStats[thisMarker].slippage+pnToSize[PnId].i1, genotypesToConsider, numOfAlleles, thisMarker.motif.size(), geomP, posSlippProb, negSlippProb);
                PnAndMarkerToGenotype[Pair<string,Marker>(PnId,thisMarker)] = changed.i1;
                Pair<bool> alleleConfidence = genotypeIsConfident(changed.i1);
                if (alleleConfidence.i1)
                    markerToStats[thisMarker].trueAlleles.insert(changed.i1.genotype.i1);
                if (alleleConfidence.i2)
                    markerToStats[thisMarker].trueAlleles.insert(changed.i1.genotype.i2);
                if (changed.i2)
                {
                    markerToStats[thisMarker].alleles.insert(changed.i1.genotype.i1);
                    markerToStats[thisMarker].alleles.insert(changed.i1.genotype.i2);
                }
                ++markerToAlleleFreqs[thisMarker].i1[changed.i1.genotype.i1];
                ++markerToAlleleFreqs[thisMarker].i1[changed.i1.genotype.i2];
                //If I am estimating the marker slippage then I should update map from Pn to labels. (before I update PnId to currentLine.PnId)
                if (!writeVcf)
                    pnToLabels[PnId] = Pair<float>(changed.i1.genotype.i1, changed.i1.genotype.i2);
                PnId = currentLine.PnId;
                clear(reads);
            }
            //Use logistic regression model to get pValue for all reads at marker
            predict_label = predict_probability(model_,probBig.x[i],prob_estimates);
            free(probBig.x[i]);
            if (i<idx)
                free(prob.x[i]);
            mapPerMarker[thisMarker][i].pValue = prob_estimates[0];
            appendValue(reads, mapPerMarker[thisMarker][i]);
        }
        free(prob.y);
        free(prob.x);
        free(probBig.x);
        free(prob_estimates);
        //Make decision about genotype for last PnId at the current marker.
        //Need to reset this set for last PN so I can add alleles from his reads
        allelesToConsider = allelesAtMarker;
        //Add alleles in reads from last PN to the set of alleles present at the marker.
        std::set<float>::iterator pnAllsEnd = PnToAlleles[PnId].end();
        for (std::set<float>::iterator pnAlls = PnToAlleles[PnId].begin(); pnAlls!=pnAllsEnd; ++pnAlls)
        {
            currAllele = *pnAlls;
            allelesToConsider.insert(currAllele);
        }
        genotypesToConsider = makeGenotypes(allelesToConsider, thisMarker.refRepeatNum).genotypes;
        changed = determineGenotype(reads, markerToStats[thisMarker].slippage+pnToSize[PnId].i1, genotypesToConsider, numOfAlleles, thisMarker.motif.size(), geomP, posSlippProb, negSlippProb);
        PnAndMarkerToGenotype[Pair<string,Marker>(PnId,thisMarker)] = changed.i1;
        Pair<bool> alleleConfidence = genotypeIsConfident(changed.i1);
        if (alleleConfidence.i1)
            markerToStats[thisMarker].trueAlleles.insert(changed.i1.genotype.i1);
        if (alleleConfidence.i2)
            markerToStats[thisMarker].trueAlleles.insert(changed.i1.genotype.i2);
        if (changed.i2)
        {
            markerToStats[thisMarker].alleles.insert(changed.i1.genotype.i1);
            markerToStats[thisMarker].alleles.insert(changed.i1.genotype.i2);
        }
        ++markerToAlleleFreqs[thisMarker].i1[changed.i1.genotype.i1];
        ++markerToAlleleFreqs[thisMarker].i1[changed.i1.genotype.i2];
        //If I am estimating the marker slippage then here is where I update the pn to labels for the last PN and write labels for this marker
        if (!writeVcf)
        {
            pnToLabels[PnId] = Pair<float>(changed.i1.genotype.i1, changed.i1.genotype.i2);
            CharString labelOutPath = options.regressionModelDirectory;
            append(labelOutPath, "/");
            append(labelOutPath, to_string(thisMarker.start));
            append(labelOutPath, "_");
            append(labelOutPath, thisMarker.motif);
            append(labelOutPath, options.iterationNumber);
            ofstream labelFile(toCString(labelOutPath));
            for (auto& label: pnToLabels)
                labelFile << label.first << "\t" << label.second.i1 << "\t" << label.second.i2 << "\n";
            pnToLabels.clear();
        }
        //Save logistic regression model to output file so I can use it in pn-slippage estimation
        CharString regressionModelDirectory = options.regressionModelDirectory;
        append(regressionModelDirectory, "/model_");
        append(regressionModelDirectory, to_string(thisMarker.start));
        append(regressionModelDirectory, "_");
        append(regressionModelDirectory, thisMarker.motif);
        const char *model_out_file = toCString(regressionModelDirectory);
        if (save_model(model_out_file,model_) != 0)
                cout << "Unable to save model for marker number " << z << endl;

        //Write vcf record for the marker I just finished.
        if (writeVcf)
        {
            //Check the number of alleles in the population, should be 2 or higher to be considered polymorphic.
            if (markerToStats[thisMarker].trueAlleles.size() < 2)
            {
                cout << "Not enough alleles at marker " << z << endl;
                mapPerMarker.erase(thisMarker);
                it = mapPerMarker.begin();
                ++z;
                continue;
            }
            //Make a String<Pair<float> > which contains a list of genotypes
            cout << "Found " << markerToStats[thisMarker].trueAlleles.size() << " alleles at marker " << z << endl;
            //Have to add ref allele to true allele set in case no one has it.
            markerToStats[thisMarker].trueAlleles.insert(thisMarker.refRepeatNum);
            genotypesAtThisMarker = makeGenotypes(markerToStats[thisMarker].trueAlleles, thisMarker.refRepeatNum);
            //First fill marker specific fields of vcfRecord
            record = fillRecordMarker(thisMarker, markerToStats[thisMarker].trueAlleles);
            //Loop over Pns and fill in PN specific fields of vcfRecord for each PN
            for (unsigned i = 0; i<length(PnIds); ++i)
            {
                thisPn = PnIds[i];
                //If a decision has been made for thisPn at thisMarker I add it to the genotypeInfos stringSet
                if ((PnAndMarkerToGenotype.count(Pair<string,Marker>(thisPn, thisMarker)) != 0))
                {
                    genotype = PnAndMarkerToGenotype[Pair<string,Marker>(thisPn, thisMarker)];
                    fillRecordPn(genotype, record, genotypesAtThisMarker,markerToStats[thisMarker].trueAlleles, thisMarker);
                    PnAndMarkerToGenotype.erase(Pair<string,Marker>(thisPn, thisMarker));
                }
                //If a decision has not been made I add a CharString with no decision(0:0,0,0,0....etc) to the set to maintain order of Pns vs genotypeInfos in output
                else
                {
                    CharString gtInfo = "./.:";
                    //Adding alleleDepth zeros
                    for (unsigned i=0; i<markerToStats[thisMarker].trueAlleles.size(); ++i)
                        append(gtInfo,"0,");
                    eraseBack(gtInfo);
                    //Adding read depth zero
                    append(gtInfo,":0:");
                    //Adding phred likelihood zeros
                    for (unsigned i=0; i<length(genotypesAtThisMarker.genotypes); ++i)
                        append(gtInfo,"0,");
                    eraseBack(gtInfo);
                    appendValue(record.genotypeInfos, gtInfo);
                }
            }
            //CharString gtInfos = concat(record.genotypeInfos);
            //cout << gtInfos << "\n" << length(record.genotypeInfos) << endl;
            ss << markerToStats[thisMarker].trueAlleles.size();
            str = ss.str();
            record.filter = ".";
            stringClear(ss,str);
            //After adding info for all PNs at thisMarker I write the record to the vcf output file
            writeRecord(out, record);
            clear(record);
        }
        markerToAlleleFreqs[thisMarker].i2 = PnsAtMarker;
        PnToAlleles.clear();
        markerSlippageOut << thisMarker.chrom << "\t" << thisMarker.start << "\t" << thisMarker.end << "\t" << thisMarker.motif << "\t" << thisMarker.refRepeatNum << "\t" << "\t" << setprecision(4) << fixed << markerToStats[thisMarker].slippage << "\t" << nAvailable << "\t" << markerToStats[thisMarker].nAlleles << "\t" << markerToStats[thisMarker].stutter << "\t" << posSlippProb << "\t" << negSlippProb << endl;
        cout << thisMarker.start << " totalSlipp: " << setprecision(4) << fixed << markerToStats[thisMarker].slippage << endl;
        cout << "Finished marker number: " << z << endl;
        mapPerMarker.erase(thisMarker);
        it = mapPerMarker.begin();
        ++z;
    }
    cout << "Finished determining genotypes" << endl;
    return 0;
}

} // namespace msGenotyper