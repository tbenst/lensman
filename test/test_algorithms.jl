using Lensman, Test
using RandomizedPropertyTest

@time @quickcheck (a+b == b+a) (a :: Int) (b :: Int)

@quickcheck (noDups(a) == noDups([r for r in eachrow(a)])) (a :: Array{Int,2})

# add 2 to group 1
badStimGroups = [
    [1  3  5  7 2],
    [2  4  6  8],
    [1  2  5  6],
    [3  4  7  8],
    [1  2  3  4],
    [5  6  7  8]
]

badGroupsPerCell = [
    [1  3  5],
    [1 2  3  5],
    [1  4  5],
    [2  4  5],
    [1  3  6],
    [2  3  6],
    [1  4  6],
    [2  4  6]
]
@test minMaxIntersection(badStimGroups, badGroupsPerCell) == (1,2)

# remove 2 from group 3
badStimGroups = [
    [1  3  5  7],
    [2  4  6  8],
    [1   5  6],
    [3  4  7  8],
    [1  2  3  4],
    [5  6  7  8]
]

badGroupsPerCell = [
    [1  3  5],
    [1 2  3  5],
    [1  4  5],
    [2  4  5],
    [1  3  6],
    [2  3  6],
    [1  4  6],
    [2  4  6]
]
@test minMaxIntersection(badStimGroups, badGroupsPerCell) == (0,1)


x = reshape(collect(1:12),4,3)
@test noDups(x)
x[1,2] = 1
@test noDups(x)==false

function validateGroupsPerCell(stimGroups, groupsPerCell, nCells, base)
    gpc = [[] for i in 1:nCells]
    for (g,group) in enumerate(stimGroups)
        for i in group
            push!(gpc[i], g)
        end
    end
    @test groupsPerCell == gpc
    nStimPerCell = Int(ceil(log(base,nCells)))
    nGroups = nStimPerCell*base
    @test length(stimGroups)==nGroups
    @test length(groupsPerCell) == nCells
    @test all(length.(groupsPerCell) .== nStimPerCell)    
end

function testSetIntersection(stimGroups, groupsPerCell)
    minIntersect, maxIntersect = minMaxIntersection(stimGroups, groupsPerCell)
    @test maxIntersect == 1
    @test minIntersect == 1
end

nCells = 8
base = 2
stimGroups, groupsPerCell = stonerStimGroups(nCells, base)
@test stimGroups[1] == collect(1:base:nCells)
@test stimGroups[2] == collect(2:base:nCells)
@test stimGroups[3] == [1,2,5,6]

validateGroupsPerCell(stimGroups, groupsPerCell, nCells, base)
testSetIntersection(stimGroups, groupsPerCell)

nCells = 64
base = 2
stimGroups, groupsPerCell = stonerStimGroups(nCells, base)
validateGroupsPerCell(stimGroups, groupsPerCell, nCells, base)
testSetIntersection(stimGroups, groupsPerCell)


nCells = 33
base = 200
stimGroups, groupsPerCell = stonerStimGroups(nCells, base)
validateGroupsPerCell(stimGroups, groupsPerCell, nCells, base)
testSetIntersection(stimGroups, groupsPerCell)

nCells = 9
base = 2
stimGroups, groupsPerCell = stonerStimGroups(nCells, base)
validateGroupsPerCell(stimGroups, groupsPerCell, nCells, base)
testSetIntersection(stimGroups, groupsPerCell)
