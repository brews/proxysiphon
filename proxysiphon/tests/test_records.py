from proxysiphon import records


def test_publication_to_citationstr():
    pub = records.Publication(authors='White, Tom; New, White',
                              published_date_or_year=1986,
                              published_title='Article title',
                              journal_name='Cool Journal',
                              volume=12, issue=3, pages=173,
                              doi='sfdjla/vcxl.3')
    goal = "White, Tom; New, White (1986): Article title. Cool Journal, 12, 3, 173, doi:sfdjla/vcxl.3"
    assert pub.to_citationstr() == goal