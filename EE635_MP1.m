% Tathagat Pal, Vikram Vasudevan

clc; clear all; close all;  % Clear the command window, workspace, and close all figures

% Parameters for Monte Carlo Simulation
snr_range = 0:1:15; % SNR values from 0 dB to 15 dB, stepping by 1 dB
num_trials = 20000;  % Number of Monte Carlo trials for each SNR value

% Define code sizes to test
code_sizes = [16, 64, 256]; % Code block sizes (N values for Polar coding)
K_factor = 0.5; % Code rate of the Polar code, K = N/2 (K is the message length)

% Initialize a matrix to store BER results for each code size and SNR value
ber_results = zeros(length(code_sizes), length(snr_range));

% Define the bit channel reliability sequence (as in 5G NR standard)
% This sequence dictates which positions in the codeword are most reliable
% The higher the index, the less reliable the bit position
Q = [0 1 2 4 8 16 32 3 5 64 9 6 17 10 18 128 12 33 65 20 256 34 24 36 7 129 66 512 11 40 68 130 ...
    19 13 48 14 72 257 21 132 35 258 26 513 80 37 25 22 136 260 264 38 514 96 67 41 144 28 69 42 ...
    516 49 74 272 160 520 288 528 192 544 70 44 131 81 50 73 15 320 133 52 23 134 384 76 137 82 56 27 ...
    97 39 259 84 138 145 261 29 43 98 515 88 140 30 146 71 262 265 161 576 45 100 640 51 148 46 75 266 273 517 104 162 ...
    53 193 152 77 164 768 268 274 518 54 83 57 521 112 135 78 289 194 85 276 522 58 168 139 99 86 60 280 89 290 529 524 ...
    196 141 101 147 176 142 530 321 31 200 90 545 292 322 532 263 149 102 105 304 296 163 92 47 267 385 546 324 208 386 150 153 ...
    165 106 55 328 536 577 548 113 154 79 269 108 578 224 166 519 552 195 270 641 523 275 580 291 59 169 560 114 277 156 87 197 ...
    116 170 61 531 525 642 281 278 526 177 293 388 91 584 769 198 172 120 201 336 62 282 143 103 178 294 93 644 202 592 323 392 ...
    297 770 107 180 151 209 284 648 94 204 298 400 608 352 325 533 155 210 305 547 300 109 184 534 537 115 167 225 326 306 772 157 ...
    656 329 110 117 212 171 776 330 226 549 538 387 308 216 416 271 279 158 337 550 672 118 332 579 540 389 173 121 553 199 784 179 ...
    228 338 312 704 390 174 554 581 393 283 122 448 353 561 203 63 340 394 527 582 556 181 295 285 232 124 205 182 643 562 286 585 ...
    299 354 211 401 185 396 344 586 645 593 535 240 206 95 327 564 800 402 356 307 301 417 213 568 832 588 186 646 404 227 896 594 ...
    418 302 649 771 360 539 111 331 214 309 188 449 217 408 609 596 551 650 229 159 420 310 541 773 610 657 333 119 600 339 218 368 ...
    652 230 391 313 450 542 334 233 555 774 175 123 658 612 341 777 220 314 424 395 673 583 355 287 183 234 125 557 660 616 342 316 ...
    241 778 563 345 452 397 403 207 674 558 785 432 357 187 236 664 624 587 780 705 126 242 565 398 346 456 358 405 303 569 244 595 ...
    189 566 676 361 706 589 215 786 647 348 419 406 464 680 801 362 590 409 570 788 597 572 219 311 708 598 601 651 421 792 802 611 ...
    602 410 231 688 653 248 369 190 364 654 659 335 480 315 221 370 613 422 425 451 614 543 235 412 343 372 775 317 222 426 453 237 ...
    559 833 804 712 834 661 808 779 617 604 433 720 816 836 347 897 243 662 454 318 675 618 898 781 376 428 665 736 567 840 625 238 ...
    359 457 399 787 591 678 434 677 349 245 458 666 620 363 127 191 782 407 436 626 571 465 681 246 707 350 599 668 790 460 249 682 ...
    573 411 803 789 709 365 440 628 689 374 423 466 793 250 371 481 574 413 603 366 468 655 900 805 615 684 710 429 794 252 373 605 ...
    848 690 713 632 482 806 427 904 414 223 663 692 835 619 472 455 796 809 714 721 837 716 864 810 606 912 722 696 377 435 817 319 ...
    621 812 484 430 838 667 488 239 378 459 622 627 437 380 818 461 496 669 679 724 841 629 351 467 438 737 251 462 442 441 469 247 ...
    683 842 738 899 670 783 849 820 728 928 791 367 901 630 685 844 633 711 253 691 824 902 686 740 850 375 444 470 483 415 485 905 ...
    795 473 634 744 852 960 865 693 797 906 715 807 474 636 694 254 717 575 913 798 811 379 697 431 607 489 866 723 486 908 718 813 ...
    476 856 839 725 698 914 752 868 819 814 439 929 490 623 671 739 916 463 843 381 497 930 821 726 961 872 492 631 729 700 443 741 ...
    845 920 382 822 851 730 498 880 742 445 471 635 932 687 903 825 500 846 745 826 732 446 962 936 475 853 867 637 907 487 695 746 ...
    828 753 854 857 504 799 255 964 909 719 477 915 638 748 944 869 491 699 754 858 478 968 383 910 815 976 870 917 727 493 873 701 ...
    931 756 860 499 731 823 922 874 918 502 933 743 760 881 494 702 921 501 876 847 992 447 733 827 934 882 937 963 747 505 855 924 ...
    734 829 965 938 884 506 749 945 966 755 859 940 830 911 871 639 888 479 946 750 969 508 861 757 970 919 875 862 758 948 977 923 ...
    972 761 877 952 495 703 935 978 883 762 503 925 878 735 993 885 939 994 980 926 764 941 967 886 831 947 507 889 984 751 942 996 ...
    971 890 509 949 973 1000 892 950 863 759 1008 510 979 953 763 974 954 879 981 982 927 995 765 956 887 985 997 986 943 891 998 766 ...
    511 988 1001 951 1002 893 975 894 1009 955 1004 1010 957 983 958 987 1012 999 1016 767 989 1003 990 1005 959 1011 1013 895 1006 1014 1017 1018 ...
    991 1020 1007 1015 1019 1021 1022 1023] + 1;

% Loop over different code sizes (N values)
for idx_N = 1:length(code_sizes)
    N = code_sizes(idx_N);  % Set the current code size (block length)
    K = K_factor * N;       % Set the message length based on code rate K/N = 0.5
    QN = Q(Q <= N);         % Select the subset of reliability sequence Q for the current block size N
    
    % For each code size, sweep across different SNR values
    for idx_snr = 1:length(snr_range)
        snr_db = snr_range(idx_snr);         % Current SNR in dB
        snr_linear = 10^(snr_db / 10);       % Convert SNR from dB to linear scale
        
        % Calculate the noise variance sigma based on the current SNR and code rate
        sigma = 1 / sqrt(2 * snr_linear * K_factor);
        
        % Initialize variables to count total bit errors and bits transmitted
        total_bit_errors = 0;
        total_bits = num_trials * K; % Total number of message bits transmitted over all trials
        
        % Monte Carlo simulation for each SNR (run num_trials times)
        for trial = 1:num_trials
            % Generate random message bits of length K
            u = randi([0 1], 1, K);
            
            % Encoder: Polar encoding using the reliability sequence Q
            d = zeros(1, N);                 % Initialize all bits to 0 (frozen bits by default)
            d(QN(N - K + 1:end)) = u;        % Assign message bits to the most reliable positions (as defined by QN)
            F = QN(1:N - K);                 % Frozen bit positions (unreliable, set to 0)
            
            
            % Polar Encoding Logic (combining the bits based on Polar transformation)
            % logic based on principles as explained in [2]
            m = 1;                           % Initialize combination factor
            for depth = log2(N) - 1:-1:0     % Go through each depth of the binary tree
                for i = 1:2*m:N              % Combine the bits pairwise
                    a = d(i:i + m - 1);
                    b = d(i + m:i + 2*m - 1);
                    d(i:i + 2*m - 1) = [mod(a + b, 2), b]; % Polar encoding logic
                end
                m = m * 2;                   % Move to the next level of the tree
            end
            codeword = d;                    % Final codeword after Polar encoding
            
            % BPSK modulation: Convert 0 -> +1, 1 -> -1
            s = 1 - 2 * codeword;
            
            % Transmit over an AWGN channel: Add Gaussian noise to the transmitted signal
            r = s + sigma * randn(1, N);     % AWGN with variance sigma^2
            
            % Successive Cancellation Decoding
            % Derived from the basic structure presented in [2]. 
            % The LLR calculations and decision rules are adjusted for our
            % specific simulation scenario. 

            L = zeros(log2(N) + 1, N);       % Log-likelihood ratio (LLR) matrix for each level
            L(1, :) = (2 * r) / sigma^2;     % Proper scaling for LLR based on received signal and noise variance
            nodestate = zeros(1, 2*N - 1);   % Node state vector to keep track of decoding decisions
            dcap = zeros(log2(N) + 1, N);    % Decoder decisions
            
            % Define the f1 and f2 functions (used for decoding)
            f1 = @(a, b) (1 - 2*(a < 0)) .* (1 - 2*(b < 0)) .* min(abs(a), abs(b)); % Min-sum function for left child
            f2 = @(a, b, c) b + (1 - 2*c) .* a; % f2 function for right child
            
            % Perform Successive Cancellation Decoding by traversing the decoding tree
            % Decoding logic from SC method, as outlined in [1], [2].
            node = 0; depth = 0; completed = 0;
            while ~completed
                if depth == log2(N) % If at the leaf nodes
                    if any(F == node + 1)
                        % Frozen bit positions are decoded as 0
                        dcap(log2(N) + 1, node + 1) = 0;
                    else
                        % Message bits: Make decision based on LLR (log-likelihood ratio)
                        dcap(log2(N) + 1, node + 1) = (L(log2(N) + 1, node + 1) < 0);
                    end
                    if node == N - 1
                        completed = 1; % If at the final node, stop decoding
                    else
                        node = floor(node / 2); depth = depth - 1;
                    end
                else
                    % Traverse through non-leaf nodes during decoding
                    npos = (2^depth - 1) + node + 1; % Position in the decoding tree
                    if nodestate(npos) == 0
                        % Go to the left child node
                        x = 2^(log2(N) - depth);
                        leftLLR = L(depth + 1, x*node + 1:x*(node + 1));
                        a = leftLLR(1:x/2); b = leftLLR(x/2 + 1:end);
                        node = node * 2; depth = depth + 1; x = x / 2;
                        L(depth + 1, x*node + 1:x*(node + 1)) = f1(a, b);
                        nodestate(npos) = 1;
                    elseif nodestate(npos) == 1
                        % Go to the right child node
                        x = 2^(log2(N) - depth);
                        leftLLR = L(depth + 1, x*node + 1:x*(node + 1));
                        a = leftLLR(1:x/2); b = leftLLR(x/2 + 1:end);
                        lnode = 2*node; ldepth = depth + 1; lx = x / 2;
                        dcapn = dcap(ldepth + 1, lx*lnode + 1:lx*(lnode + 1));
                        node = node * 2 + 1; depth = depth + 1; x = x / 2;
                        L(depth + 1, x*node + 1:x*(node + 1)) = f2(a, b, dcapn);
                        nodestate(npos) = 2;
                    else
                        % Pass information up to the parent node
                        x = 2^(log2(N) - depth);
                        lnode = 2 * node; rnode = 2 * node + 1;
                        cdepth = depth + 1; cx = x / 2;
                        dcapl = dcap(cdepth + 1, cx*lnode + 1:cx*(lnode + 1));
                        dcapr = dcap(cdepth + 1, cx*rnode + 1:cx*(rnode + 1));
                        dcap(depth + 1, x*node + 1:x*(node + 1)) = [mod(dcapl + dcapr, 2), dcapr];
                        node = floor(node / 2); depth = depth - 1;
                    end
                end
            end
            
            % Extract the decoded message from the reliable positions
            decoded_msg = dcap(log2(N) + 1, QN(N - K + 1:end));
            
            % Calculate the number of bit errors in this trial by comparing with transmitted message
            bit_errors_in_trial = sum(decoded_msg ~= u);
            total_bit_errors = total_bit_errors + bit_errors_in_trial;
        end
        
        % Compute the Bit Error Rate (BER) for the current SNR and code size
        ber_results(idx_N, idx_snr) = total_bit_errors / total_bits;
    end
end

% Plot the BER vs. SNR curve for each code size
figure;
for idx_N = 1:length(code_sizes)
    semilogy(snr_range, ber_results(idx_N, :), 'DisplayName', sprintf('N = %d', code_sizes(idx_N)), 'LineWidth', 2);
    hold on;
end

xlabel('SNR (dB)', 'FontSize', 14);
ylabel('Bit Error Rate (BER)', 'FontSize', 14);
title('Polar Codes BER vs SNR', 'FontSize', 14);
legend('show', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 14);

ber_results

% The structure of this code is based on material from several sources.
% Specifically, the basic polar encoding and decoding logic were inspired by 
% code examples from the following resources:
% References:
% [1] https://github.com/vetrvl/polar-codes-list/tree/master
% [2] https://www.youtube.com/watch?v=rB0rhQKyV34&list=PLyqSpQzTE6M81HJ26ZaNv0V3ROBrcv-Kc&index=27


% The BER decreases as the SNR increases, which is expected. Higher SNR 
% values indicate better signal quality, resulting in fewer bit errors 
% during transmission. This downward trend is consistent across 
% different block sizes (N = 16, N = 64, N = 256).

% For low SNR values (0 to 2 dB), smaller block sizes (N = 16) perform better 
% than larger block sizes (N = 64 and N = 256). This happens because larger code 
% sizes tend to have more complex encoding/decoding processes, which require 
% a higher SNR to fully benefit from their error-correcting capabilities. 
% At low SNR, larger code sizes do not provide significant gains in error 
% correction and instead introduce additional overhead.

% The plot stops after 5 dB likely because, at higher SNR values, 
% the simulations resulted in an extremely low or zero BER (no bit errors), 
% making it difficult to gather meaningful data for SNR values beyond 5 dB. 
% At these high SNRs, the communication channel becomes very clean, and 
% polar codes perform well. No errors are observed within simulation's limited trials.
