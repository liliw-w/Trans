input=$1
while IFS= read -r line
do
  echo "$line"
  bash -c $line
done < "$input"
