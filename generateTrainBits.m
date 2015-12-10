trainBits50 = randi([0 1], 1, 50);
trainBits75 = randi([0 1], 1, 75);
trainBits100 = randi([0 1], 1, 100);
save('trainBits', 'trainBits50', 'trainBits100', 'trainBits150');