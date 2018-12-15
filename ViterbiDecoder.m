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
  returnValue=(y1-x1)^2+(y2-x2)^2;
end

function returnValue = ber(sourceA,sourceB)

if length(sourceA) ~= length(sourceB)
    error -1;
end

returnValue = sum(sourceA ~= sourceB);
end

%This encoder is from lab3 and does not always generate a "correct" signal for the viterbi decoder
function returnValue = convolutionalEncoder(source)
  extendedSource = [0 0 source];
  returnValue = zeros (1, length(source)*2);
  j=1;
  for i=3:1:length(extendedSource)
    returnValue(j)=mod(extendedSource(i-2)+extendedSource(i), 2);
    returnValue(j+1)=mod(extendedSource(i-2)+extendedSource(i-1)+extendedSource(i), 2);
    j=j+2; 
  end
end

function returnValueMatrix = viterbiLikelyhood(returnValueMatrix,value,source,currentColumnPosition,currentPosition,outputValue,possibleSteps,finish,HD)
if(HD=='HARD')
switchValues =  [0 0, 0 1, 1 0,1 1];
else
switchValues =  mapping([0 0, 0 1, 1 0,1 1]);
endif
sdValue=0;
switch (value)  
  %00
  case [switchValues(1) switchValues(2)]
    sdValue=returnValueMatrix(1,currentColumnPosition-1)+SD(source(currentPosition),outputValue(1,1),source(currentPosition+1),outputValue(1,2));
    if(returnValueMatrix(possibleSteps(1,1),currentColumnPosition)>sdValue)
      returnValueMatrix(possibleSteps(1,1),currentColumnPosition)=sdValue;
    endif
    if (finish==false)
      sdValue=returnValueMatrix(1,currentColumnPosition-1)+SD(source(currentPosition),outputValue(1,3),source(currentPosition+1),outputValue(1,4));
      if(returnValueMatrix(possibleSteps(1,2),currentColumnPosition)>sdValue)
        returnValueMatrix(possibleSteps(1,2),currentColumnPosition)=sdValue;
      endif
    endif
    %01
  case [switchValues(3) switchValues(4)]
    sdValue=returnValueMatrix(2,currentColumnPosition-1)+SD(source(currentPosition),outputValue(2,1),source(currentPosition+1),outputValue(2,2));
    if(returnValueMatrix(possibleSteps(2,1),currentColumnPosition)>sdValue)
      returnValueMatrix(possibleSteps(2,1),currentColumnPosition)=sdValue;
    endif
    if (finish == false)
      sdValue=returnValueMatrix(2,currentColumnPosition-1)+SD(source(currentPosition),outputValue(2,3),source(currentPosition+1),outputValue(2,4));
      if(returnValueMatrix(possibleSteps(2,2),currentColumnPosition)>sdValue)
      returnValueMatrix(possibleSteps(2,2),currentColumnPosition)=sdValue;
      endif
    endif
  %10
  case [switchValues(5) switchValues(6)]
    sdValue=returnValueMatrix(3,currentColumnPosition-1)+SD(source(currentPosition),outputValue(3,1),source(currentPosition+1),outputValue(3,2));
    if(returnValueMatrix(possibleSteps(3,1),currentColumnPosition)>sdValue)
      returnValueMatrix(possibleSteps(3,1),currentColumnPosition)=sdValue;
    endif
    if (finish == false)
      sdValue=returnValueMatrix(3,currentColumnPosition-1)+SD(source(currentPosition),outputValue(3,3),source(currentPosition+1),outputValue(3,4)) ; 
      if(returnValueMatrix(possibleSteps(3,2),currentColumnPosition)>sdValue)
        returnValueMatrix(possibleSteps(3,2),currentColumnPosition)=sdValue;
      endif
    endif
  %11
  case [switchValues(7) switchValues(8)]
    sdValue=returnValueMatrix(4,currentColumnPosition-1)+SD(source(currentPosition),outputValue(4,1),source(currentPosition+1),outputValue(4,2));
    if(returnValueMatrix(possibleSteps(4,1),currentColumnPosition)>sdValue)
      returnValueMatrix(possibleSteps(4,1),currentColumnPosition)=sdValue;
    endif
    if (finish==false)
      sdValue=returnValueMatrix(4,currentColumnPosition-1)+SD(source(currentPosition),outputValue(4,3),source(currentPosition+1),outputValue(4,4));
      if(returnValueMatrix(possibleSteps(4,2),currentColumnPosition)>sdValue)
        returnValueMatrix(possibleSteps(4,2),currentColumnPosition)=sdValue;
      endif
    endif
  otherwise 
    error(-1);
endswitch
end

function returnValue = getValues(position,HD)
  switch (position)
    case 1
      if(HD=='HARD')
        returnValue=[0 0];
      else
        returnValue=[1 1];
      endif      
    case 2
      if(HD=='HARD')
        returnValue=[0 1];
      else
        returnValue=[1 -1];
      endif
    case 3
       if(HD=='HARD')
        returnValue=[1 0];
      else
       returnValue=[-1 1];
      endif
    case 4
       if(HD=='HARD')
        returnValue=[1 1];
      else
        returnValue=[-1 -1];
      endif
  endswitch
end

function returnValue = getBitValue(currentRowPosition,previousRowPosition)
  switch (previousRowPosition)
    case 1
      if(currentRowPosition==1)
        returnValue = 0;
      else
        returnValue = 1;
      endif
    case 2
      if(currentRowPosition==3)
        returnValue = 0;
      else
        returnValue = 1;
      endif
    case 3
      if(currentRowPosition==1)
        returnValue = 0;
      else
        returnValue = 1;
      endif
    case 4
      if(currentRowPosition==3)
        returnValue = 0;
      else
        returnValue = 1;
      endif
    otherwise
      error(-1)
    endswitch
end


function [returnValue,stepMatrix,diff] = mostLiklyStep(stepMatrix,returnValueMatrix,startPoint,possibleSteps,diff=0); 
    currentColumn=returnValueMatrix(:,startPoint);
    rowsWithSmallesValue=find(currentColumn==min(currentColumn));
    if(length(rowsWithSmallesValue)==1)
      if(startPoint==1)
        stepMatrix(1,1)=0;
        diff=0;
        returnValue=getBitValue(rowsWithSmallesValue,1);
      else
        if(diff!=0)
          previousRowsWithSmallesValue=find(returnValueMatrix(:,startPoint-diff)==min(returnValueMatrix(:,startPoint-diff)));
          for i=1:length(previousRowsWithSmallesValue)
            if(any (possibleSteps(previousRowsWithSmallesValue(i),:) == rowsWithSmallesValue))
              stepMatrix(previousRowsWithSmallesValue(i),startPoint-diff)=0; 
              returnValue=getBitValue(previousRowsWithSmallesValue(i),find(stepMatrix(:,startPoint-(1+diff))==0));
              return;
            endif
          end
        endif
        stepMatrix(rowsWithSmallesValue,startPoint-diff)=0 ;      
        returnValue=getBitValue(rowsWithSmallesValue,find(stepMatrix(:,startPoint-1)==0));
        diff=0;
      endif
    else
      diff=diff+1;
      [returnValue,stepMatrixReturn]=mostLiklyStep(stepMatrix,returnValueMatrix,startPoint+1,possibleSteps,diff);
       stepMatrix=stepMatrixReturn;
    endif  
end

function returnValue = mostLiklySteps(len,returnValueMatrix,possibleSteps)
  returnValues=zeros(1,len-3);
  len = length(returnValueMatrix);
  stepMatrix=ones(4,len-2)*9999;
  for i=1:1:len-2
    [returnValues(i),stepMatrixReturn]=mostLiklyStep(stepMatrix,returnValueMatrix,i,possibleSteps);
    stepMatrix=stepMatrixReturn;
end 
returnValue=returnValues;
end

function returnValue = viterbiDecoderFiveSeven(source,HD='HARD')
  len=length(source)/2;
  possibleSteps=[1 2;3 4;1 2;3 4];
  extendedSource=[0 0 source 0 0];
  outputValue=[0 0 1 1;0 1 1 0;1 1 0 0;1 0 0 1];
  returnValue=zeros(len+2);
  extendedSource=[0 0 source 0 0 0 0];
  if (HD!='HARD')
    outputValue=mapping(outputValue);
  endif
  returnValueMatrix=ones(4,len+3)*9999;
  currentSteps=possibleSteps(1,:);
  returnValueMatrix(1,1)=0;
  currentPosition=1;
  currentColumnPosition=1;
  for k=1:len
    steps=find(returnValueMatrix(:,currentColumnPosition)<9999);
    currentPosition=currentPosition+2;
    currentColumnPosition=currentColumnPosition+1;
    for i=1:length(steps)
      returnValueMatrix=viterbiLikelyhood(returnValueMatrix,getValues(i,HD),extendedSource,currentColumnPosition,currentPosition,outputValue,possibleSteps,false,HD);
    end
  end
  
    %finish
    for k=0:1:1
      steps=find(returnValueMatrix(:,currentColumnPosition)<9999);
      currentPosition=currentPosition+2;
      currentColumnPosition=currentColumnPosition+1;
      for i=1:length(steps)
        returnValueMatrix=viterbiLikelyhood(returnValueMatrix,getValues(i,HD),extendedSource,currentColumnPosition,currentPosition,outputValue,possibleSteps,true,HD);
      end
  end
  returnValue = mostLiklySteps(len,returnValueMatrix,possibleSteps)(2:len+1); 
end


infoBits = round(rand(1,10000));
%infoBits=[1 1 1 0 1 0 0 1]
%% see comment at convolutinalEncoder
convolutionalEncodedSignal=mapping(convolutionalEncoder(infoBits));
%y=[1 1 0 0 1 0 0 1 0 0]
viterbiDecodedSignal=viterbiDecoderFiveSeven(convolutionalEncodedSignal,'SOFT');
BER = ber(infoBits,viterbiDecodedSignal)





