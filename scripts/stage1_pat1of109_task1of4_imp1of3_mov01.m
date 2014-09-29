clc;close all;clear all;

load('d:\Dropbox\Signals\EEG Motor Movement-Imagery Dataset\Processed\S001\R03\S001R03')

T = floor(sqrt(min(size(mov1{1},2))));
F3 = 32;

signal = mov1{F3}(1,:);
port = perPort(signal,T);

