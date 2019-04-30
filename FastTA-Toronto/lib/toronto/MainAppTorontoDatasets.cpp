#include <eo>

#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "testset/TestSetDescription.h"
#include "init/ETTPInit.h"

#include "neighbourhood/ETTPneighborEval.h"
#include "eval/eoETTPEval.h"
#include "neighbourhood/ETTPneighbor.h"
#include "toronto/TorontoTestSet.h"
#include "algorithms/eo/eoCellularEARing.h"
#include "algorithms/eo/eoCellularEAMatrix.h"
#include "algorithms/eo/Mutation.h"
#include "algorithms/eo/Crossover.h"

//#include "eoSimpleGA.h"
//#include "algorithms/eo/eoGenerationContinue.h"
#include "algorithms/eo/eoGenerationContinuePopVector.h"
#include "eoSelectOne.h"
#include "algorithms/eo/eoSelectBestOne.h"
//#include <eoDetSelect.h>
#include "algorithms/eo/eoDeterministicTournamentSelectorPointer.h" // eoDeterministicTournamentSelector using boost::shared_ptr

// For counting the # evaluations
#include "eval/eoNumberEvalsCounter.h"
#include "eval/eoETTPEvalNumberEvalsCounter.h"
#include "eval/statistics/eoETTPEvalWithStatistics.h"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "statistics/ExamMoveStatistics.h"
#include "statistics/optimised/ExamMoveStatisticsOpt.h"



using namespace std;



#define MAINAPP_DEBUG


extern int getSANumberEvaluations(double tmax, double r, double k, double tmin);
//
// Rong Qu's exam evaluator
//
extern void runExamEvaluator(string const& _torontoDatasetPath, string const& _outputDir, string const& _torontoDatasetName);



// These function is defined below
void runCellularEA(string const& _outputDir, TestSet const& _testSet);
void generateExamMoveStatistics(const string &_outputDir, const TestSet &_testSet);


void runTA(TestSet const& _testSet, string const& _outputDir,
           moSimpleCoolingSchedule<eoChromosome> &_coolSchedule);



void runTorontoDatasets(int _datasetIndex, string const& _testBenchmarksDir, string const& _outputDir) {

    vector<TestSetDescription> torontoTestSet;
    //
    // To evaluate the algorithms, we use Version I of Toronto benchmarks
    // TestSetDescription info: <name>, <description>, <# time slots>
    //
    torontoTestSet.push_back(TestSetDescription("car-f-92", "Carleton University",                   32));
    torontoTestSet.push_back(TestSetDescription("car-s-91", "Carleton University",                   35));
    torontoTestSet.push_back(TestSetDescription("ear-f-83", "Earl Haig Collegiate",                  24));
    torontoTestSet.push_back(TestSetDescription("hec-s-92", "Ecole des Hautes Etudes Commerciales",  18));
    torontoTestSet.push_back(TestSetDescription("kfu-s-93", "King Fahd University",                  20));
    torontoTestSet.push_back(TestSetDescription("lse-f-91", "London School of Economics",            18));
    torontoTestSet.push_back(TestSetDescription("pur-s-93", "Purdue University",                     42));
    torontoTestSet.push_back(TestSetDescription("rye-s-93", "Ryerson University",                    23));
    torontoTestSet.push_back(TestSetDescription("sta-f-83", "St. Andrews High school",               13));
    torontoTestSet.push_back(TestSetDescription("tre-s-92", "Trent University",                      23));
    torontoTestSet.push_back(TestSetDescription("uta-s-92", "University of Toronto, Arts & Science", 35));
    torontoTestSet.push_back(TestSetDescription("ute-s-92", "University of Toronto, Engineering",    10));
    torontoTestSet.push_back(TestSetDescription("yor-f-83", "York Mills Collegiate",                 21));
    /////////////////////////////////////////////////////////////////////////////////////////////////
//    copy(torontoTestSet.begin(), torontoTestSet.end(), ostream_iterator<TestSetDescription>(cout, "\n"));

    vector<TestSetDescription>::iterator it = torontoTestSet.begin() + _datasetIndex;
    // Create TestSet instance
    TorontoTestSet testSet((*it).getName(), (*it).getDescription(), _testBenchmarksDir);
    // Load dataset
    TorontoTestSet* ptr = &testSet;
    ptr->load();
    // Fixed length timetables
    int numPeriods = (*it).getPeriods();
    // Set the # periods in the timetable problem data
    ptr->getTimetableProblemData()->setNumPeriods(numPeriods);
    // Print timetable problem data
//    cout << *(ptr->getTimetableProblemData().get()) << endl;
    // Print testset info
//    cout << testSet << endl;

    // Generate exam move statistics
    generateExamMoveStatistics(_outputDir, testSet);

    // Run test set
//    runCellularEA(_outputDir, testSet);


    // Parameters: string const& _torontoDatasetPath, string const& _outputDir, string const& _torontoDatasetName
//    runExamEvaluator(testSet.getRootDirectory()+"/", _outputDir, testSet.getName());

    ///////////
    /// TA
    ///
//    moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.1, 0.001, 5, 2e-5);//  Uta 3.68 1h // LIGHT COOL SCHEDULE
//    moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.1, 0,000001, 5, 2e-5); // INTENSIVE COOL SCHEDULE
//    runTA(testSet, _outputDir, coolSchedule);

    ///////////


}


// Generate exam move statistics
void generateExamMoveStatistics(string const& _outputDir, TestSet const& _testSet) {
    ///////////////////////////////////////////////////
    //
    // Generate exam move statistics
    //
    ///////////////////////////////////////////////////
    cout << "///////////////////////////////////////////////////" << endl;
    cout << "//" << endl;
    cout << "// Generate exam move statistics" << endl;
    cout << "//" << endl;
    cout << "///////////////////////////////////////////////////" << endl;
    // Number of thresholds
    int numBins = 10;
    // TA cooling schedule
    //    moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.1, 0.99, 3, 2e-5); // Sch #0
    //      moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.5, 0.001, 5, 2e-5); // Sch #11

//    moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.1, 0.001, 5, 2e-5); // Light Cooling schedule - USED

//          moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.1, 0.0001, 5, 2e-5); // Sch #2, 6m30s

//        moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.1, 0.00001, 5, 2e-5); // Uta 3.13, Sch #3 - FastTA Pur (4.48) 30 min
        //    moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.5, 0.00001, 5, 2e-5); // Uta 3..., Sch #4

//        moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.1, 0.000001, 5, 2e-5); // Uta 3.03, Sch #5

    // Mid-way cooling schedule
    moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.1, 0.00001, 5, 2e-5); // USED

    // Intensive cooling schedule
//    moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.1, 0.0000005, 5, 2e-5); // USED

    //    moSimpleCoolingSchedule<eoChromosome> coolSchedule(0.1, 0.0000001, 5, 2e-5); // Sch #6

    ExamMoveStatistics examMoveStats(_testSet, _outputDir, numBins, coolSchedule);
    examMoveStats.run(_testSet, _outputDir);

//    ExamMoveStatisticsOpt examMoveStats(_testSet, _outputDir, numBins, coolSchedule);
//    examMoveStats.run();

    /////////////////////////////////////////////////////
    ///
    /// NOTE:
    /// For creating the output file containing threshold information and
    /// Creating the output file containing information of move counts in each threshold per exam
    /// USE ExamMoveStatistics instead of ExamMoveStatistics*Opt*
    ///
    /////////////////////////////////////////////////////

    cout << "END" << endl;
}



void runTA(TestSet const& _testSet, string const& _outputDir,
           moSimpleCoolingSchedule<eoChromosome> &_coolSchedule) {

    // Creating the output filename
    stringstream sstream;
    sstream << _outputDir << "/TA_" << _testSet.getName() << "_cool_"
            << _coolSchedule.initT << "_" << _coolSchedule.alpha << "_"
            << _coolSchedule.span << "_" << _coolSchedule.finalT << ".txt";
    string outFilename;
    sstream >> outFilename;
    std::cout << outFilename << std::endl;
    // Output file
    std::ofstream outFile(outFilename);
    ///////////////////////////////////////////////////////////
    long maxNumEval;
    double tmax = _coolSchedule.initT, r = _coolSchedule.alpha,
           k = _coolSchedule.span, tmin = _coolSchedule.finalT;
    // max # evaluations
    maxNumEval = getSANumberEvaluations(tmax, r, k, tmin);
    std::cout << "numberEvaluations = " << maxNumEval << std::endl;
    // Print max # evaluations to file
    outFile << "numberEvaluations = " << maxNumEval << std::endl;
    ///////////////////////////////////////////////////////////
    // Solution initializer
    ETTPInit<eoChromosome> init(*_testSet.getTimetableProblemData().get());
    // Generate initial solution
    eoChromosome initialSolution;
    init(initialSolution);
    // # evaluations counter
    eoNumberEvalsCounter numEvalsCounter;
    // eoETTPEval used to evaluate the solutions; receives as argument an
    // eoNumberEvalsCounter for counting neigbour # evaluations
    eoETTPEvalNumberEvalsCounter<eoChromosome> fullEval(numEvalsCounter);
    // Evaluate solution
    fullEval(initialSolution);

    //
    // Local search used: Threshold Accepting algorithm
    //
    // moTA parameters
    boost::shared_ptr<ETTPKempeChainHeuristic<eoChromosome> > kempeChainHeuristic(
                new ETTPKempeChainHeuristic<eoChromosome>());
    ETTPneighborhood<eoChromosome>neighborhood(kempeChainHeuristic);
    // ETTPneighborEval which receives as argument an
    // eoNumberEvalsCounter for counting neigbour # evaluations
    ETTPneighborEvalNumEvalsCounter<eoChromosome> neighEval(numEvalsCounter);

    moTA<ETTPneighbor<eoChromosome> > ta(neighborhood, fullEval, neighEval, _coolSchedule);

    /////// Write to output File ///////////////////////////////////////////
    cout << "Start Date/Time = " << currentDateTime() << endl;
    // Write Start time and algorithm parameters to file
    outFile << "Start Date/Time = " << currentDateTime() << endl;
    outFile << "TA parameters:" << endl;
    outFile << "cooling schedule: " << _coolSchedule.initT << ", " << _coolSchedule.alpha << ", "
            << _coolSchedule.span << ", " << _coolSchedule.finalT << endl;
    outFile << _testSet << std::endl;

    /////////////////////////////////////////
    // Get current time
    time_t now;
    double seconds = 0.0;
    time(&now);  // get current time; same as: now = time(NULL)
    /////////////////////////////////////////

    cout << "Before TA - initialSolution.fitness() = " << initialSolution.fitness() << endl;

    // Apply TA to the solution
    ta(initialSolution);

    // Validate solution
//    initialSolution.validate();

    cout << "After TA - initialSolution.fitness() = " << initialSolution.fitness() << endl;

    // Write best solution to file
    outFile << "Solution fitness = " << initialSolution.fitness() << endl;
    // Print real # evaluations performed
    std::cout << "# evaluations performed = " << numEvalsCounter.getTotalNumEvals() << std::endl;
    outFile << "# evaluations performed = " << numEvalsCounter.getTotalNumEvals() << endl;
    // Print solution timetable to file
    outFile << initialSolution << endl;

    //
    // Create solution file and run Rong Qu's exam evaluator
    //
    string solutionFilename = _testSet.getName() + ".sol";
    ofstream solutionFile(_outputDir + solutionFilename);
    // Write best solution to file
    initialSolution.printToFile(solutionFile);
    // Close solution file
    solutionFile.close();
    // Parameters: string const& _torontoDatasetPath, string const& _outputDir, string const& _torontoDatasetName
    runExamEvaluator(_testSet.getRootDirectory()+"/", _outputDir, _testSet.getName());

    /////////////////////////////////////////
    // Get current time
    time_t final;
    time(&final);  // get current time
    // Get difference in seconds
    seconds = difftime(final, now);
    /////////////////////////////////////////
    cout << "End Date/Time = " << currentDateTime() << endl;
    cout << "Seconds elapsed = " << seconds << endl;
    // Write to file
    outFile << "End Date/Time = " << currentDateTime() << endl;
    outFile << "Seconds elapsed = " << seconds << endl;
}









