function [vpcode] = getVpCode
% getVpCode prompts for researcher initials and participant number
% and compiles a 4-character participant code (e.g., 'AA01').

FlushEvents('keyDown');

% Ask for researcher initials
researcherInitials = input('\n\n>>>> Please type researcher initials: ', 's');
researcherInitials = upper(strtrim(researcherInitials)); % Convert to uppercase & trim spaces

% Ask for participant number
participantNumber = input('\n>>>> Please type researcher participant number: ', 's');
participantNumber = strtrim(participantNumber); % Trim spaces

% Ensure participant number is two digits
if length(participantNumber) == 1
    participantNumber = strcat('0', participantNumber);
end

% Construct vpcode
vpcode = strcat(researcherInitials, participantNumber);

end