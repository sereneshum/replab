function message = lawsMessage(errorDesc, context, assertArgs, testVarNames, testVarValues)
    if ~isequal(context, '')
        message = [errorDesc char(10) context];
    else
        message = errorDesc;
    end
    for i = 1:length(assertArgs)
        try
            assertArgs{i} = str(assertArgs{i});
        catch
            assertArgs{i} = moxunit_util_elem2str(assertArgs{i});
        end
    end
    message = sprintf(message, assertArgs{:});
    message = [message char(10) 'with variables' newline];
    maxLength = max(cellfun(@(x) length(x), testVarNames));
    for i = 1:length(testVarNames)
        vn = testVarNames{i};
        vn = [blanks(maxLength - length(vn)) vn ':' char(10)];
        try
            vv = str(testVarValues{i});
        catch
            vv = moxunit_util_elem2str(testVarValues{i});
        end
        message = [message sprintf('%s\n%s\n\n', vn, vv)];
    end
end
