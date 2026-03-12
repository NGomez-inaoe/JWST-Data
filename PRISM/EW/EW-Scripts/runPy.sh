for f in Ew-Scripts/*.py; do
    python "$f" &
done
wait
