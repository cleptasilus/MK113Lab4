clc
clear all

% mapping 0 to 1V and 1 to -1V 
function returnValue = mapping(source) 
  source(source>0)=-1;
  source(source==0)=1;
  returnValue = source;
end

% demapping the signal
function returnValue = demapping(source)
  source(source>=0)=0;
  source(source<0)=1;
  returnValue=source;
end

function returnValue = SD(y1,x1,y2,x2)
  [y1 y2]
  [x1 x2]
  returnValue=(y1-x1)^2+(y2-x2)^2
end

function returnValueMatrix = viterbiLikelyhood(returnValueMatrix,value,source,currentColumnPosition,currentPosition,outputValue,possibleSteps,finish)
%switchValues =  mapping([0 0, 0 1, 1 0,1 1])
switchValues =  [0 0, 0 1, 1 0,1 1];
value
sdValue=0;
switch (value)  
  %00
  case [switchValues(1) switchValues(2)]
    sdValue=returnValueMatrix(1,currentColumnPosition-1)+SD(source(currentPosition),outputValue(1,1),source(currentPosition+1),outputValue(1,2))
    if(returnValueMatrix(possibleSteps(1,1),currentColumnPosition)>sdValue)
      returnValueMatrix(possibleSteps(1,1),currentColumnPosition)=sdValue
    endif
    if (finish==false)
      sdValue=returnValueMatrix(1,currentColumnPosition-1)+SD(source(currentPosition),outputValue(1,3),source(currentPosition+1),outputValue(1,4))
      if(returnValueMatrix(possibleSteps(1,2),currentColumnPosition)>sdValue)
        returnValueMatrix(possibleSteps(1,2),currentColumnPosition)=sdValue
      endif
    endif
    %01
  case [switchValues(3) switchValues(4)]
    sdValue=returnValueMatrix(2,currentColumnPosition-1)+SD(source(currentPosition),outputValue(2,1),source(currentPosition+1),outputValue(2,2))
    if(returnValueMatrix(possibleSteps(2,1),currentColumnPosition)>sdValue)
      returnValueMatrix(possibleSteps(2,1),currentColumnPosition)=sdValue
    endif
    if (finish == false)
      sdValue=returnValueMatrix(2,currentColumnPosition-1)+SD(source(currentPosition),outputValue(2,3),source(currentPosition+1),outputValue(2,4))
      if(returnValueMatrix(possibleSteps(2,2),currentColumnPosition)>sdValue)
      returnValueMatrix(possibleSteps(2,2),currentColumnPosition)=sdValue
      endif
    endif
  %10
  case [switchValues(5) switchValues(6)]
    sdValue=returnValueMatrix(3,currentColumnPosition-1)+SD(source(currentPosition),outputValue(3,1),source(currentPosition+1),outputValue(3,2))
    if(returnValueMatrix(possibleSteps(3,1),currentColumnPosition)>sdValue)
      returnValueMatrix(possibleSteps(3,1),currentColumnPosition)=sdValue
    endif
    if (finish == false)
      sdValue=returnValueMatrix(3,currentColumnPosition-1)+SD(source(currentPosition),outputValue(3,3),source(currentPosition+1),outputValue(3,4))  
      if(returnValueMatrix(possibleSteps(3,2),currentColumnPosition)>sdValue)
        returnValueMatrix(possibleSteps(3,2),currentColumnPosition)=sdValue
      endif
    endif
  %11
  case [switchValues(7) switchValues(8)]
    sdValue=returnValueMatrix(4,currentColumnPosition-1)+SD(source(currentPosition),outputValue(4,1),source(currentPosition+1),outputValue(4,2))
    if(returnValueMatrix(possibleSteps(4,1),currentColumnPosition)>sdValue)
      returnValueMatrix(possibleSteps(4,1),currentColumnPosition)=sdValue
    endif
    if (finish==false)
      sdValue=returnValueMatrix(4,currentColumnPosition-1)+SD(source(currentPosition),outputValue(4,3),source(currentPosition+1),outputValue(4,4))
      if(returnValueMatrix(possibleSteps(4,2),currentColumnPosition)>sdValue)
        returnValueMatrix(possibleSteps(4,2),currentColumnPosition)=sdValue
      endif
    endif
  otherwise 
    error(-1);
endswitch
end

function returnValue = getValues(position)
  switch (position)
    case 1
%      returnValue=[1 1]
      returnValue=[0 0];
    case 2
%      returnValue=[1 -1]
      returnValue=[0 1];
    case 3
%      returnValue=[-1 1]
      returnValue=[1 0];
    case 4
%      returnValue=[-1 -1]
      returnValue=[1 1];
  endswitch
end

function returnValue = viterbiDecoderFiveSeven(source)
  
possibleSteps=[1 2;3 4;1 2;3 4]

%outputValue=mapping([0 0 1 1;0 1 1 0;0 1 1 0;0 0 1 1])
%source=mapping([1 1,0 1,1 0,0 1,0 1,0 0,0 0])
outputValue=[0 0 1 1;0 1 1 0;0 1 1 0;0 0 1 1]
%source=[1 1,0 1,1 0,0 1,0 1,0 0,0 0]
  returnValue=zeros(7);
  extendedSource=[0 0 source];
  returnValueMatrix=ones(4,7)*9999;
  currentSteps=possibleSteps(1,:)
  returnValueMatrix(1,1)=0;
  currentPosition=1;
  currentColumnPosition=1;
  for k=1:5
    steps=find(returnValueMatrix(:,currentColumnPosition)<9999)
    currentPosition=currentPosition+2;
    currentColumnPosition=currentColumnPosition+1;
    for i=1:length(steps)
      returnValueMatrix=viterbiLikelyhood(returnValueMatrix,getValues(i),extendedSource,currentColumnPosition,currentPosition,outputValue,possibleSteps,false)
    end
  end
  
  %finish for now hard coded to test
      steps=find(returnValueMatrix(:,currentColumnPosition)<9999)
    currentPosition=currentPosition+2;
    currentColumnPosition=currentColumnPosition+1;
    for i=1:length(steps)
      returnValueMatrix=viterbiLikelyhood(returnValueMatrix,getValues(i),extendedSource,currentColumnPosition,currentPosition,outputValue,possibleSteps,true)
    end
  
  
end


source=[1 1,0 0,1 0,0 1,0 0,0 0]
source=[1 1,1 0,1 0,0 1,0 0,0 0]
viterbiDecoderFiveSeven(source);





