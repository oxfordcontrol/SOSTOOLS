function display(var)
% display(var) prints the dpvar object var using the polynomial class
% format. This function is automatically called when a dpvar output
% statement is not terminated by semicolon.
% 
% INPUTS:
% var: dpvar
% 
% OUTPUTS:
% command line display of elements of dpvar
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ, MP, SS - 07/30/2021
% Replace "sum(.,'all')" with "sum(sum(.))" for pre-R2018b, DJ - 02/07/2022


% tol indicates the max number of command-line rows that the output is
% allowed to occupy without user confirmation
% Set to 0 for default value (width of command window)
tol = 0;

% Output will be displayed under name varname
if isempty(inputname(1))
    varname = 'ans';
else
    varname = inputname(1);
end

if isempty(var)
% % % If the object is empty, inform the user of this
    if ~isequal(get(0,'FormatSpacing'),'compact')
        disp(' ');
    end
    disp([varname ' =']);
    disp([' Empty dpvar: ',num2str(size(var,1)),'-by-',num2str(size(var,2))]);
    if ~isequal(get(0,'FormatSpacing'),'compact')
        disp(' ');
    end

elseif var.chkval   
% % % A valid (nonempty) dpvar gets displayed in polynomial format % % %
    
    % % Setup: convert all matrix elements to chars (which we can display)
    
    % Convert to polynomial
    P = dpvar2poly(var);
    
    % Collect elements of polynomial matrix as chars in a cell
    Pchar = char(P);
    szP = size(P);
    
    % Check the size required to display each element of the matrix...
    nchar = cellfun('size',Pchar,2);
    % ... and the max size required for each column of the matrix...
    maxr = max(nchar,[],1);
    % ... to determine the width in the command window required to display
    % the object in full
    numc = sum(maxr);

    
    % % Display the actual object
    % Introduce a vspace if required
    if ~isequal(get(0,'FormatSpacing'),'compact')
        disp(' ');
    end
   
    rootprops = fieldnames(get(0));
    % Establish the maximum number of chars that can fit on one row in the command window
    if any(strcmp(rootprops,'CommandWindowSize'))
        ws = get(0,'CommandWindowSize');
    else
        % Older versions of matlab don't allow access to the
        % window size.  Modify window size (ws) if polynomial
        % line breaks are not in the right spot.
        ws = 80;
    end
    
    % Account for spacing between matrix columns, and start and end spacing
    maxchar = max([1, ws(1)-6-3*(szP(2)-1)]); 
    if tol==0
        tol = ws;  % Max number of rows to print in command window
    end
    
    % Display the object as full matrix if possible, or elementwise if not
    if numc<=maxchar
        % Each row of the object can fit in one line of the command window
        if szP(1)<=tol
            % The object is not too big --> we can savely display
            poly_display_rows(varname,Pchar,szP,nchar);
        else
            % The object is (very) big --> check if the user wishes to continue
            msg = ['Warning: Your object will require a large number of rows'...
                   ' (',num2str(sum(szP(1))),') to display in your command window.'];
            fprintf(1,[msg,'\n'])
            userinp = input('Do you wish to continue with display? Y or N?\n','s');
            userinp = strrep(userinp,' ','');
            if strcmpi(userinp,'N') || strcmpi(userinp,'No')
                return
            end
            poly_display_rows(varname,Pchar,szP,nchar);
        end
    else
        % Each element of the object gets a separate row
        nrows = ceil(nchar./maxchar);
        if sum(sum(nrows+2))<=tol
            % The object is not too big --> we can savely display
            poly_display_elems(varname,Pchar,szP,ws);
        else
            % The object is (very) big --> check if the user wishes to continue
            msg = ['Warning: Your object will require a large number of rows'...
                   ' (estimated ',num2str(sum(sum(nrows+2))),') to display in your command window.'];
            fprintf(1,[msg,'\n'])
            userinp = input('Do you wish to continue with display? Y or N?\n','s');
            userinp = strrep(userinp,' ','');
            if strcmpi(userinp,'N') || strcmpi(userinp,'No')
                return
            end
            poly_display_elems(varname,Pchar,szP,ws);
        end
    end
    
    % Introduce another vspace if required
    if ~isequal(get(0,'FormatSpacing'),'compact')
        disp(' ');
    end
    
else % var.chkval = 0   
    error('Your object is not a valid dpvar!')
    
end

end



function poly_display_rows(name,val,szP,nchar)
% Displays polynomial elements of the cell "val", under the name "name".
% Polynomial is displayed in full matrix format.
% "szP" must be the size of the associated matrix-valued polynomial, and
% "nchar" an integer cell array specfiying the number of characters 
% required to display each element of "val".
% Copied from @polynomial display, by PJS

maxr = max(nchar,[],1);
disp([name ' = ']);
for i1 = 1:szP(1)
    if all(szP==[1 1])
        d = '  ';
    else
        d = '  [ ';
    end
    for i2 = 1:szP(2)
        d = [d blanks(maxr(i2)-nchar(i1,i2)) val{i1,i2}];
        if i2~=szP(2)
            d = [d ', '];
        elseif ~all(szP==[1 1])
            d = [d ']'];
        end
    end
    disp(d);
end

end



function poly_display_elems(name,val,szP,ws)
% Displays polynomial elements of the cell "val", under the name "name".
% Polynomial is displayed entry-by-entry going down columns.
% "szP" must be the size of the associated matrix-valued polynomial, and
% "ws" an integer specifying the allowed number of characters in a single
% row of the command window.
% Copied from @polynomial display, by PJS

for i2 = 1:szP(2) %i1 is the row and i2 is the column
    if min(szP) > 1
        disp(['----  Column ' int2str(i2) ' ----------'])
        disp(' ')
    end
    for i1 = 1:szP(1)
        % Display the (i1,i2) entry
        if ~all(szP==[1 1])
            disp([name '(' int2str(i1) ',' int2str(i2) ')  = ']);
        else
            disp([name ' = ']);
        end
        
        % Break lines at +,-, or *
        sr = val{i1,i2};
        while ~isempty(sr) %length(sr)>0
            if length(sr) < (ws(1)-6)
                disp(['  ' sr]);
                sr = [];
            else
                idx1 = sort([strfind(sr,'-') strfind(sr,'+') strfind(sr,'*')]);
                idx1 = [setdiff(idx1,1) length(sr)];
                %idx2 = max(find(idx1< (ws(1)-6) ));
                idx2 = find(idx1< (ws(1)-6), 1, 'last' );
                if isempty(idx2)
                    disp(['  ' sr(1:idx1(1)-1)]);
                    sr = sr(idx1(1):end);
                else
                    disp(['  ' sr(1:idx1(idx2)-1)]);
                    sr = sr(idx1(idx2):end);
                end
            end
        end
        if ~all([i1 i2]==szP)
            disp(' ');
        end
    end
end

end