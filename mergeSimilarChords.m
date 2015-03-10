% this is a process to merge chords with same bass and similar type into
% one chord. Assume there are two chords to be merged, the mergee is the one with
% less notes, while the merger is the one with more notes; if both are with
% the same number of notes, the second one is to be merged to the first one
% there are totally three super chord types: 1. those with 3 (major
% type), 2. those with 3b (minor type), 3. those without 3 or 3b (universal
% type). Chords within same type can be merged. Type 1 and Type 2 can be
% merged into type 3 or merge from type 3, but type 1 and type 2 cannot
% be merged into each other.
function [chordogram, bassgram, treblegram, boundaries] = mergeSimilarChords(chordogram, bassgram, treblegram, boundaries, chordmode)



