__author__ = 'chwheele'
import numpy, random
from processingFunctions import merge_dicts


def createInitialRandPopulation(rangeDict, initPopNum, intFlags=None):
    """
    Create a random initial population based on the values in range Dict and the desired number of
    population members (sets of parameters to test).
    :param rangeDict (dict): The keys of this dictionary are the parameter names being tested.
    Each key corresponds to a tuple that has the minimum and maximum values allowed for this parameter in the
    search (minParam, maxParam).
    :param initPopNum (int): The number of population members to return (sets of parameters to test.)
    :param intFlags (list of strings): These are the keys in params Dict that will be flagged to have integer values.
    The default value type is a float.
    :return initPopulation (list of Dicts): This returns a list of dictionaries, each corresponding to a single 
    population member. Each dictionary will have the same "keys" as the rangeDict, while each value in the dictionary 
    will be between the numbers specified in the range dict.
    """
    # this list will hold the population dictionaries. 
    initPopulation = []
    # int flags should be a list of strings. These strings as the keys in the 
    # rangeDict that must have integer values. In the case where None a specified we make this an empty set.
    if intFlags is None:
        intFlags = []
    # loop for every member of the population
    for popIndex in range(initPopNum):
        # make a fresh dictionary for every population member
        popMemberDict = {}
        # loop through the parameter keys of the rangeDict
        for key in rangeDict.keys():
            (minValue, maxValue) = rangeDict[key]
            # make sure the maxValue is bigger than or equal to the min value
            if maxValue < minValue:
                tempVar = minValue
                minValue = maxValue
                maxValue = tempVar
                del tempVar
            # calculate the min max difference
            diff = maxValue - minValue
            if maxValue == minValue:
                #case where the maxValue and minValue are thee same so there is only one possible parameter
                popMemberDict[key] = minValue
            else:
                # Case where minValue < maxValue
                # If statement Below: True if this value should a integer, False if it should be a float
                if key in intFlags:
                    # saves value as an integer
                    popMemberDict[key] = int(numpy.round((random.random() * diff) + minValue))
                else:
                    # saves value as a float
                    popMemberDict[key] = (random.random() * diff) + minValue
        # save the new population member (set of parameters to test) after exiting the parameter loop. 
        initPopulation.append(popMemberDict)
    return initPopulation


def rotate(popMatrix, seedIndexArray, rotationAmount):
    """
    :param popMatrix (narray): an order matrix of population member's parameters
    :param seedIndexArray (narray): an array of uniquie indexes for an array of len popNum
    :param rotationAmount: The number of indexes of which the seed array (and thus the popMatrix) will be rotated.
    :return rotPopMembers (narray): a array of population members randomized by the seed array, and then rotated,
    by the integer stored in rotation amount.
    """
    popNum = len(seedIndexArray)
    # add a rotation amound to the seed indexes
    rotatedIndexs = seedIndexArray + rotationAmount
    # trim the indexes that are greater than the number of population members
    modulusRotatedIndexes = numpy.mod(rotatedIndexs, popNum)
    rotPopMembers = popMatrix[modulusRotatedIndexes, :]
    return rotPopMembers


def differentialEvolution(parentPop, rangeDict, crossOverProb=0.5, strategy=2, F=0.8):
    """
    This makes a new generation from an existing generation using differential evolutions

    :param parentPop (list of dictionaries): Should be a list of dictionaries order from best performances to worst,
        where each dictionary corresponding to a single population member.
        Each dictionary will have the same "keys" as the rangeDict, while each value in the
        dictionary will be between the numbers specified in the range dict.
    :param rangeDict (list of dictionaries): The keys of this dictionary are the parameter names being tested.
        Each key corresponds to a tuple that has the minimum and maximum values allowed for this parameter in the
        search (minParam, maxParam).
    :param crossOverProb (0 <= float < 1): probability of genes (parameter values), crossing over from the parent
        population and replacing the mutant gene. Set to 0 to not allow crossover, setting to 1 means
        all genes will be crossovers from the parents. Setting to 1 is not desired and effectively clones the parent
        population into the next generation.
    :param strategy (int in (1, 2, 3, 4, 5)): This are possible differential evolution strategies.
        Strategies 1 and 4 use the best population member are the starting point for all mutations,
        These give quick convergence to a single solotion
    :param F (0 <= float): Practically this should be less than 1 to avoid large leaps. This is a scalar that
        multiplies the difference vector between the genes of 2 or more parents. The closer to 0, the smaller the
        effect of gene mutation (smaller changes in parameters per step).
    :return: 
    """
    ### This gets vector math intensive, so initialize some variables if you need this code to run faster.

    ## Initialization not important for python
    # popMemberList1=numpy.zeros((popNum,parameterNum))
    # popMemberList2=numpy.zeros((popNum,parameterNum))
    # popMemberList3=numpy.zeros((popNum,parameterNum))
    # popMemberList4=numpy.zeros((popNum,parameterNum))
    # popMemberList5=numpy.zeros((popNum,parameterNum))
    # intermediateVectors=numpy.zeros((popNum,parameterNum))
    # maskIntermediateVectors=numpy.zeros((popNum,parameterNum), dtype=bool)
    # maskIntermediateVectors2=numpy.zeros((parameterNum,popNum))
    # maskIntermediateVectors3=numpy.zeros((parameterNum,popNum))

    ## Initialization not important for python
    # maskIntermediateVectors4=numpy.zeros((popNum,parameterNum))
    # mask4oldPop=numpy.zeros((popNum,parameterNum))

    ## Initialization not important for python
    # someOtherRotators = numpy.zeros(popNum)
    # rotatingArray4ExponentialCrossover=numpy.zeros(parameterNum)
    # indexArray1=numpy.zeros(popNum)
    # indexArray2=numpy.zeros(popNum)
    # indexArray3=numpy.zeros(popNum)
    # indexArray4=numpy.zeros(popNum)
    # indexArray5=numpy.zeros(popNum)
    # indexPointerArray=numpy.zeros(4)
    # iterationNum=1

    """
    This is the start of the parameter unwrapping sequence. This mostly house keeping to convert from convient
    dictionaries to arrays for faster (an more intuitive for Caleb) math.  
    """
    mutantPopMatrix = None
    # these are used to create various arrays
    popNum = len(parentPop)
    parameterNum = len(rangeDict)
    # convert the parent population list of dictionaries a matrix for numpy
    parameterKeys = rangeDict.keys()
    popMatrix = numpy.zeros((popNum,parameterNum))
    for (popIndex, popMember) in list(enumerate(parentPop)):
        for (paramIndex, parameter) in list(enumerate(parameterKeys)):
            try:
                popMatrix[popIndex, paramIndex] = popMember[parameter]
            except:
                singleTonList = popMember[parameter]
                popMatrix[popIndex, paramIndex] = singleTonList[0]

    """
    Convert the best member (expected to be the first value of parentPop)
    to a population matrix of only the best member (is is good to be king).
    """
    popOfOnlyBestMember=numpy.zeros((popNum,parameterNum))
    bestMember = parentPop[0]
    for popIndex in range(popNum):
        for (paramIndex, parameter) in list(enumerate(parameterKeys)):
            try:
               popOfOnlyBestMember[popIndex, paramIndex] = bestMember[parameter]
            except:
                singleTonList=bestMember[parameter]
                popOfOnlyBestMember[popIndex, paramIndex] = singleTonList[0]



    """
    Here we make some arrays to shuffle the population matrix around.
    These shuffled arrays switch around the 'breeding partners' of the population members 
    (set of parameters to test).
    """
    # This makes sure the rotation is at least 1 and the possible rotations are all unique
    # This is to ensure that population goes un-rotated
    popRotationIndexArray = numpy.arange(1, popNum)
    random.shuffle(popRotationIndexArray)

    """
    seedIndexArray contains a randomized array of unique indexes corresponding to each different population members.
    The values in popMatrix are expected to be ordered from best to worst, by suffling the seed array we make use it
    to a randomized list of population members.
    """
    seedIndexArray = numpy.arange(popNum)
    random.shuffle(seedIndexArray)


    """
    This makes 'n = numberOfRandPopulationRotations' lists of the population members.
    Each list is randomized compared to the ordered population matrix.
    Each list is a unique rotation of elements of all lists. This ensures that for a given parameter in each list
    will be from a unique parent. This solves the problem of accidental asexual reproduction within the ordinal
    version of this code
    """
    numberOfRandPopulationRotations = 5
    # this dictionary stores the rotated lists of the population members.
    randomPopOrderDict = {}
    for n in range(numberOfRandPopulationRotations):
        # each step adds a unique rotated list to the dictionary for use in the strategy section of the code.
        randomPopOrderDict['popMemberList' + str(n)] \
            = rotate(popMatrix, seedIndexArray, rotationAmount=popRotationIndexArray[n], popNum=popNum)


    ### Now we will make the cross-over matrix
    # Make a mask that is True randomly at a frequency given by crossOverProb
    # True means no crossover for this gene, mutant genes only
    maskIntermediateVectors = crossOverProb < numpy.random.rand(popNum, parameterNum)
    # True means parent gene makes a crossover, into the mutant population
    # this makes the inverse mask to maskIntermediateVector
    mask4oldPop = numpy.logical_not(maskIntermediateVectors)

    """
    Here are examples of 5 possible means differential evolution (5 DE stratagies)
    """
    # Strategy 1
    if (strategy == 1):                  # DE/best/1
        mutantPopMatrix = popOfOnlyBestMember + F * (randomPopOrderDict['popMemberList1'] - randomPopOrderDict['popMemberList2']) # differential variation
        mutantPopMatrix = randomPopOrderDict['popMemberList0'] * mask4oldPop + mutantPopMatrix * maskIntermediateVectors # crossover
    # Strategy 2
    elif (strategy == 2):                  # DE/rand/1
        mutantPopMatrix = randomPopOrderDict['popMemberList3'] + F * (randomPopOrderDict['popMemberList1'] - randomPopOrderDict['popMemberList2']) # differential variation
        mutantPopMatrix = randomPopOrderDict['popMemberList0'] * mask4oldPop + mutantPopMatrix*maskIntermediateVectors # crossover
    # Strategy 3
    elif (strategy == 3):                  # DE/rand-to-best/1
        mutantPopMatrix = randomPopOrderDict['popMemberList0'] + F * (popOfOnlyBestMember - randomPopOrderDict['popMemberList0']) + F * (randomPopOrderDict['popMemberList1'] - randomPopOrderDict['popMemberList2'])
        mutantPopMatrix = randomPopOrderDict['popMemberList0'] * mask4oldPop + mutantPopMatrix*maskIntermediateVectors # crossover
    # Strategy 4
    elif (strategy == 4):                  # DE/best/2
        mutantPopMatrix = popOfOnlyBestMember + F * (randomPopOrderDict['popMemberList1'] - randomPopOrderDict['popMemberList2'] + randomPopOrderDict['popMemberList3'] - randomPopOrderDict['popMemberList4']) # differential variation
        mutantPopMatrix = randomPopOrderDict['popMemberList0'] * mask4oldPop + mutantPopMatrix * maskIntermediateVectors # crossover
    # Strategy 5
    elif (strategy == 5):                  # DE/rand/2
        mutantPopMatrix = randomPopOrderDict['popMemberList5'] + F * (randomPopOrderDict['popMemberList1'] - randomPopOrderDict['popMemberList2'] + randomPopOrderDict['popMemberList3'] - randomPopOrderDict['popMemberList4']) # differential variation
        mutantPopMatrix = randomPopOrderDict['popMemberList0'] * mask4oldPop + mutantPopMatrix * maskIntermediateVectors # crossover

    """
    Here we do some housekeeping to make sure all the mutations fall with the ranges specified in the rangeDict
    """
    # get a min and a max for each parameter.
    parameterMinList=[]
    parameterMaxList=[]
    for (paramIndex,parameter) in list(enumerate(parameterKeys)):
        (parameterMin,parameterMax) = rangeDict[parameter]
        parameterMinList.append(min([parameterMin,parameterMax]))
        parameterMaxList.append(max([parameterMin,parameterMax]))

    mutantPop = []
    for popIndex in range(popNum):
        popMember = {}
        memberParams=mutantPopMatrix[popIndex,:]
        for (paramIndex,parameter) in list(enumerate(parameterKeys)):
            # make sure the parameters are in the range specified by 'rangeDict'
            parameter2test = memberParams[paramIndex]
            parameterMin = parameterMinList[paramIndex]
            parameterMax = parameterMaxList[paramIndex]
            # force it to be within the minimum
            while parameter2test <= parameterMin:
                diff = abs(parameter2test < parameterMin)
                parameter2test = parameterMin + (1.5 * diff)
            # force it to be with in the maximum,
            while parameterMax <= parameter2test:
                diff = abs(parameter2test-parameterMax)
                parameter2test = parameterMax - (1.5 * diff)
            popMember[parameter]= parameter2test
        mutantPop.append(popMember)

    """
    Now more housekeeping to change things back into conveinte dictionaries for other parts of the code.
    """
    cleanedParentPop=[] # this version of the population has gotten rid of singleton lists
    for popIndex in range(popNum):
        popMember = {}
        memberParams=popMatrix[popIndex,:]
        for (paramIndex,parameter) in list(enumerate(parameterKeys)):
            popMember[parameter]=memberParams[paramIndex]
        cleanedParentPop.append(popMember)

    combinedParentPop = []
    for popIndex in range(popNum):
        oldParentPopMember = parentPop[popIndex]
        cleanedParentPopMember = cleanedParentPop[popIndex]

        combinedParentPop.append(merge_dicts(oldParentPopMember, cleanedParentPopMember))


    return mutantPop, combinedParentPop


if __name__ == "__main__":
    """
    This is a test to demonstrate the above code. It only runs when diffEvo.py is run directly as in 
    'python diffEvo.py' but does not when called from another python file, for example, importing a definition.
    """
    numberOfPopulationMembers = 50
    numberOfParametersPerPopMember = 5 # used to make rangeDict
    crossOverProb = 0.1
    strategy = 1
    F = 0.8

    # make a test range dictionary, to specify the bounds of each parameter.
    rangeDict = {"testVal" + str(n): (0.0 + (10.0 * float(n)), 1.0 + (10.0 * float(n))) for n in range(numberOfParametersPerPopMember)}

    # make the list of population member dictionaries that are made of sets of parameters.
    parentPop = createInitialRandPopulation(rangeDict,
                                            initPopNum=numberOfPopulationMembers,
                                            intFlags=None)
    
    # Do a single step of differential evolution
    mutantPop, combinedParentPop = differentialEvolution(parentPop,
                                                         rangeDict,
                                                         crossOverProb=crossOverProb,
                                                         strategy=strategy,
                                                         F=F)
    print "The parent population is:\n", combinedParentPop
    print "The mutant population is:\n", mutantPop
