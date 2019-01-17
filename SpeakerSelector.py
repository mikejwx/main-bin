import random
# first number is the number of days since last spoke
# second number is 1 if staff, 2 if student
my_list = ['Kirsty Hanley']           * (766 * 1)**2 \
        + ['Carol Halliwell']         * (766 * 1)**2 \
        + ['Peter Clark']             * (738 * 1)**2 \
        + ['Oscar Martinez-Alvarado'] * (731 * 1)**2 \
        + ['Lewis Blunn']             * (38  * 2)**2 \
        + ['Nigel Roberts']           * (654 * 1)**2 \
        + ['Miguel Teixeira']         * (640 * 1)**2 \
        + ['Jacob Maddison']          * (10  * 2)**2 \
        + ['Bob Plant']               * (598 * 1)**2 \
        + ['Ben Harvey']              * (563 * 1)**2 \
        + ['Thorwald Stein']          * (556 * 1)**2 \
        + ['Julia Curio']             * (542 * 1)**2 \
        + ['Benjamin Courtier']       * (528 * 2)**2 \
        + ['Alex Baker']              * (521 * 1)**2 \
        + ['Bethan Harris']           * (10  * 2)**2 \
        + ['Humphrey Lean']           * (381 * 1)**2 \
        + ['Kaja Milczewska']         * (3   * 2)**2 \
        + ['James Gilmore']           * (353 * 2)**2 \
        + ['Ambrogio Volonte']        * (346 * 1)**2 \
        + ['Matthew Priestley']       * (339 * 2)**2 \
        + ['Alec Vessey']             * (332 * 2)**2 \
        + ['Dan Shipley']             * (297 * 2)**2 \
        + ['Will McIntyre']           * (24  * 2)**2 \
        + ['Ryan Williams']           * (241 * 2)**2 \
        + ['Mark Muezelfeldt']        * (0   * 2)**2 \
        + ['Michael Johnston']        * (3   * 2)**2 \
        + ['Helen Dacre']             * (185 * 1)**2 \
        + ['Chris Holloway']          * (178 * 1)**2 \
        + ['Sue Gray']                * (171 * 1)**2

i = 0
chosen = []
while i < 14:
    new_speaker = random.choice(my_list)
    if new_speaker not in chosen:
        print str(i + 1) + '. ' + new_speaker + '\n'
        chosen.append(new_speaker)
        i += 1