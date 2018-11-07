for j in `find /Users/rjspaargaren/Documents/Thesis/Code_Dan/spider-dev/Thesis_data/lid/* -name "n*" -type d -exec basename {} \;`;
do echo $j;
folname=$j;
cd /Users/rjspaargaren/Documents/Thesis/Code_Dan/spider-dev/Thesis_data/lid/$folname;
mkdir output;
spider -options_file bu_input.opts;
done;
