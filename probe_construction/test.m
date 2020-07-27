clear;clc;
regionProps = zeros(6,8);
regionProps(1,:) = [1,5,27,35,57,78,82,104];
regionProps(2,:) = 20;
padLength = 0;
selectedRegionData = TRDesigner.TileSplitRegions(regionProps, padLength);
expected = [5,27,82,104];
fprintf('Expected v.s. returned:\n')
disp(expected)
disp(selectedRegionData(1,:))
