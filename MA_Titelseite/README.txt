Usage

1. Save the file  "MA_Titelseite.sty" in the directory in which the master's thesis is saved.
2. Add the following data to the header of your thesis and change the values of the variables.


%Name of the author of the thesis 
\authornew{X Y}
%Date of birth of the Author
\geburtsdatum{ }
%Place of Birth
\geburtsort{New York, U.S.A.}
%Date of submission of the thesis
\date{ }

%Name of the Advisor
% z.B.: Prof. Dr. Peter Koepke
\betreuer{Advisor: Prof. Dr. X Y}
%name of the second advisor of the thesis
\zweitgutachter{Second Advisor: Prof. Dr. X Y}

%Name of the Insitute of the advisor
%z.B.: Mathematisches Institut
\institut{Mathematisches Institut}
%\institut{Institut f\"ur Angewandte Mathematik}
%\institut{Institut f\"ur Numerische Simulation}
%\institut{Forschungsinstitut f\"ur Diskrete Mathematik}
%Title of the thesis 
\title{This is only an example}
%Do not change!
\ausarbeitungstyp{Master's Thesis  Mathematics}

3. With "\maketitle"  you can now generate the title page.
